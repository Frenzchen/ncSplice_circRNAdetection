#!/usr/bin/env ruby

module Alignment

	extend self
	
	# Trim end of alignment if mismatch is in the last 2 bp, likely to be a mapping error.
	#
	#
	# array	- Array of Integers from mapping with breakpoint_downstream or
	#					breakpoint_upstream methods.
	# mm		- Number of allowed mismatches.
	#
	# Return integer of trimmed alignment array.
	def trim(array, mm)
		array[-8..-1].count(0) == 0 ? array.length : 0
		array = array[0..-3] if array[-2..-1] == [0,1] && mm == 1
		array = array[0..-2] if array[-1] == 0 && mm == 1
		array[-8..-1].count(0) == 0 ? array.length : 0
	end


	# Extend 3' anchor to upstream genomic region.
	# Use the trim-method to remove mismatches from the alignment end.
	#
	#
	# ref - Reference sequence to which DNA sequence will be compared.
	# mm  - Number of allowed mismatches.
	#
	# Return integer of alignment length.
	def upstream(seq, reference, mm)
		a = []
		mismatch = 0

		(seq.length - 1).downto(0).each do |i|
			if seq[i] == reference[i]
				a << 1
			else
				mismatch += 1
				return trim(a, mm) if mismatch > mm
				a << 0
			end
		end
	end


	# Extend 5' anchor to downstream genomic region.
	# Use the trim-method to remove mismatches from the alignment end.
	#
	#
	# ref - Reference sequence to which DNA sequence will be compared.
	# mm  - Number of allowed mismatches.
	#
	# Return integer of alignment length.
	def downstream(seq, reference, mm)
		a = []
		mismatch = 0
	
		0.upto(seq.length - 1).each do |i|
			if seq[i] == reference[i]
				a << 1
			else
				mismatch += 1
				return trim(a, mm) if mismatch > mm
				a << 0
			end
		end
	end


	# Create reverse complement of DNA.
	#
	#
	# dna - DNA sequence
	#
	# Return reverse complement.
	def reverse_complement(dna)
		complement = []
		dna.each_char do |s|
			if s == 'A'
				complement << 'T'
			elsif s == 'T'
				complement << 'A'
			elsif s == 'G'
				complement << 'C'
			elsif s == 'C'
				complement << 'G'
			elsif s == 'N'
				complement << 'N'
			end	
		end
		complement.join('').reverse
	end
	
	
	# Check read quality based on phred value.
	#
	#
	# quality	- Obserevd phred value of read.
	# phred		- Minimal phred value allowed.
	#
	# Return boolean.
	def quality_ok?(quality, phred)
  	quality.each_char.all? {|x| (x.ord - 33) >= phred}
	end


	# Get read orientation.
	# Orientation can't be defined if library type is unstranded.
	# If library-type is set to fr-firststrand, R1 determines the strand.
	# If library-type is set to fr-secondstrand, R2 determines the strand.
	#
	#
	# library_type	- Library type ('fr-unstranded', 'fr-firstrand' or 'fr-secondstrand').
	# mate					- Read mate (R1 or R2)
	# strand				- Strand to which mate was mapped (1 or -1).
	#
	# Returns orientation of read as 1 or -1 (+ strand or - strand).
	def read_orientation(library_type, mate, strand)

		if library_type == 'fr-unstranded'
			orientation = strand
		else 
			# fr-firststrand
			if mate == '1'
				strand == 1 ? orientation = -1 : orientation = 1	
			elsif mate == '2'
				strand == 1 ? orientation = 1 : orientation = -1	
			end
			
			# fr-secondstrand
			orientation = (orientation * -1) if library_type == 'fr-secondstrand'
			orientation
		end
	end
	
	
	# Compare number of mismatches reported by MD:Z tag to allowed number of mismatches.
	#
	#
	# mdz - String of MD:Z tag from bowtie2 alignment.
	# mm  - Integer of max. number of allowed mismatches.
	#
	# Return boolean.
	def max_mismatches?(mdz, mm)
		mdz = mdz.split(':').last
		counter = 0
		
		if !mdz.nil?
			%w[A T G C].each {|x| counter += mdz.count(x) }
			counter <= mm
		else
			FALSE
		end
	end
	
	
	# Report genompic distance over which read mapped.
	#
	#
	# cigar				- Cigar string from bowtie2 alignment.
	# read_length	- Length of sequecning reas.
	#
	# Return false if insertions etc present, otherwise genomic length as interger.
	def genomic_mappinglength(cigar, read_length)
		if cigar == "#{read_length}M"
			read_length
		elsif cigar.include?('N') && %w(I D S H P).all? {|x| !cigar.include?(x)}
			cigar = cigar.split(/N|M/).collect {|x| x.to_i}
			cigar.inject {|sum, x| sum + x}
		else
			false
		end
	end
	
	
	# Check whether spliced read has mapped mate in the right position and orientation.
	#
	#
	# mate1	- Mapping information for mate 1.
	# mate2	- Mapping information for mate 2.
	#
	# Return false if insertions etc present, otherwise genomic length as interger.
	def paired?(mate1, mate2)
		mate1_chr, mate1_start, mate1_stop, mate1_strand = mate1
		mate2_chr, mate2_start, mate2_stop, mate2_strand = mate2
		conditions = []
		
		conditions.push(mate1_chr == mate2_chr)
		conditions.push(mate1_start.to_i <= mate2_start && mate1_stop.to_i >= mate2_stop)
		conditions.push(mate1_strand.to_i != mate2_strand)
		conditions.all? { |con| con == true }
	end
	
	
	# Generate possible kmers from overhanging DNA.
	# First and last basepair are reported separately, the rest with 2bp sliding window.
	#
	#
	# dna	- Overhanging DNA.
	#
	# Returns array.
	def kmers(dna)
		motifs = []
		0.upto(dna.length - 2) { |i| motifs << dna[i..i+1] }
		motifs
	end  	

	
	# Compare motif via 2bp sliding window to know U2 and U12 splice sites.
	# Each motif pair (upstream + downstream) is scored, based on whether it is canonical or
	# 	whether it is 1 or 2 base pairs longer ("stronger" if you have two basepairs).
	# Canonical & 2 base pairs: + 1
	# Canonical & 1 base pairs: + 0.75
	# Non-canonical & 2 base pairs: + 0.5
	# Non-canonical & 1 base pairs: + 0.25
	#
	#
	# overhang_dna	- Overhanging dna (left/right of read extendedable to both sides) 
	# orientation		- Orientation of strand (1: sense, -1, anti-sense)
	#
	# Returns integer corrsponding to index of highest scored motif
	def score_motifs(up_dna, down_dna, orientation)

		overhang = up_dna.length

		if orientation == -1
			up = Alignment.reverse_complement(down_dna)
			down = Alignment.reverse_complement(up_dna)
		else
			up = up_dna
			down = down_dna
		end

  	motifs = {:don => ['GT', 'GC', 'AT', 'CT'], :acc => ['AG', 'AC']}
  	score = {:don => Array.new(overhang+1, 0), :acc => Array.new(overhang+1, 0)}
		
		uKmers = kmers(up)
		dKmers = kmers(down)
		
		# scores, different depending on whether it is canonical and length 
		# higher score (more "power") if you have two base pairs to play
		# downstream (donor) motif
		dKmers.each_with_index do |k, i|
			if motifs[:don].include?(k)
				k == 'GT' ? score[:don][i] += 1.0 : score[:don][i] += 0.5
			end
		end
  	
  	# upstream (acceptor) motif	
		uKmers.each_with_index do |k, i|
			if motifs[:acc].include?(k)
				k == 'AG' ? score[:acc][i] += 1.0 : score[:acc][i] += 0.5
			end
		end

  	# find maximal motif score	
  	summed_scores = score.values.transpose.map {|x| x.reduce(:+)}
  	max_score = summed_scores.max 
		several_max = summed_scores.count(max_score)
		several = summed_scores.select {|x| x > 0}.length
		
  	if max_score > 0

  		i = summed_scores.index {|x| x == max_score}

 			if motifs[:acc].include?(uKmers[i]) && motifs[:don].include?(dKmers[i])
 				motif = dKmers[i] + ':' + uKmers[i]
 			elsif motifs[:acc].include?(uKmers[i]) && !motifs[:don].include?(dKmers[i])
 				motif =  '-:' + uKmers[i]	
 			elsif !motifs[:acc].include?(uKmers[i]) && motifs[:don].include?(dKmers[i])
 				motif = dKmers[i] + ':-'
 			end
  		
  		{:index => i, :scores => [score[:acc].join('|'), score[:don].join('|')],\
  		 :no_scores => [several_max, several].join('|'), :max_score => max_score, \
  		 :motif => motif}

  		# unknown motif	
  	else
  		{:index => '-', :scores => [score[:acc].join('|'), score[:don].join('|')],\
  		 :no_scores => '0|0', :max_score => 0, :motif => '-'}
  	end
	end
	
	
	# Correct breakpoints after motif trimming.
	#
	# 
	# motif_summary	-	Motif summary hash.
	# upstream_bp		- Original upstream breakpoint.
	# downstream_bp	- Original downstream breakpoint.
	# overhang			- Overhang of read.
	#
	# upstream correction
	# upstream_bp + (index - 2) + 2
	# (index - 2): 2 additional nt upstream to define splice motif -> shift in index
	# +2 nt afterwards to get to first exon position
	# so +2/-2 cancel each other out
	#
	# downstream correction	
	# downstream_bp - overlap + (index - 2) + 2
	# - overlap: get nt from wich we have to add the remaining nt of "overlap" 	
	# (index - 2): 2 additional nt upstream to define splice motif -> shift in index
	# +2 nt afterwards to get to last exon position
	# so +2/-2 cancel each other out
	#
	# Returns corrected breakpoints (integer hash).
	def correct_breakpoint(motif_summary, upstream_bp, downstream_bp, overhang, strand)
		if motif_summary[:index].kind_of?(Integer)
			i = motif_summary[:index]
			if strand == 1
				upstream_bp = upstream_bp + i
				downstream_bp = downstream_bp - overhang + i
			else
				upstream_bp = upstream_bp + overhang - i
				downstream_bp = downstream_bp - i
			end
		end
		[upstream_bp, downstream_bp]
	end
end
