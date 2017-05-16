#!/usr/bin/env ruby

# 13/03/2016
# version: ruby 2.0.0
# 
# Wrapper script for circRNA detection.
#
#
# Options (update!):
# u - Fastq-file with unmapped reads from tophat in bam-format.
# p - Prefix used for all output files.
# x - Path to bowtie2-index directory and base name <index_directory>/<bt2-index>.
# a - Length of anchors, 20 bp is recommended.
# l - Length of read.
# f - Path to directory with one fasta files for each chromosome.
# s - File with chromosomes that should excluded.
#			One chromosome per line.
# singletones - File with paired reads for which only one mate mapped.
# sequencing-type - SE or PE for single-end or paired-end.
#
# Remarks:
# Implement additional option to delete intermediate files.


require 'optparse'
require 'open3'
require_relative "../lib/alignments.rb"
require_relative "../lib/analysis.rb"
require_relative "../lib/bamClass.rb"

# options
##########################################################################################

options = {}
optparse = OptionParser.new do |opts|
  opts.banner = "ncSplice version 0.0.0 by Franziska Gruhl (franziska.gruhl@swiss.ch)\n
  Usage: ncSplice_prepareUnmapped.rb -u <unmapped reads> -f <file format: bam, fasta or fastq> -p <prefix> --sequencing-type <se or pe> -x <bowtie-index>/<index-prefix> [further options]"

  opts.on('-h', '--help', 'Display help screen.') do
    puts opts
    exit
  end
	
	opts.on('-v', '--version', 'Print ncSplice version and dependencies.') do
		puts '# ncSplice'
		puts 'ncSplice v0.1.0'
		puts "\n# Dependencies"
		puts ['ruby >= 2.0.0', 'samtools >= 1.0.0', 'bowtie2 >= 2.1.0'].join("\n")
    exit
  end
	
	options[:phred] = 25
	options[:anchor] = 20
	options[:bowtie_version] = 'bowtie'
	
	opts.on('-u', '--input <filename>', String, 'An input file with unmapped reads. The file format (bam, fastq or fasta) needs to be specified via the -f option.') {|u| options[:unmapped] = u}
	opts.on('-f', '--format <string>', String, 'The input format of unmapped reads file, can be fastq, bam or fasta.') {|f| options[:format] = f}
	opts.on('-q', '--quality <integer>', Integer, 'The minimal phred quality each base in the unmapped read needs to have to be considered for further analysis. The default is set to 25.') {|q| options[:phred] = q}
  opts.on('-p', '--prefix <string>', String, 'A prefix or label, which will be used for all output files.') {|p| options[:prefix] = p}  
  opts.on('--sequencing-type <string>', String, 'The sequencing type. Use \'se\' for single-end and \'pe\' for paired-end sequencing.') {|seq| options[:sequencing] = seq}
  opts.on('-a', '--anchor-length <integer>', Integer, 'The length of the read anchor used in the remapping procedue, the default is 20 bp. Shorter anchors will decrease the mapping precision and longer anchors will cause a reduction in the number of candidates.') {|a| options[:anchor] = a}
  opts.on('-b', '--bowtie-version <String>', String, 'The bowtie version to use, can be either bowtie or bowtie2 (default bowtie).') {|b| options[:bowtie_version] = b}
  opts.on('-x', '--bowtie-index <directory>', String, 'The bowtie-index diretory and base name: <index-directory>/<bt2-index>.') {|x| options[:bowtie_index] = x}
  opts.on('-m', '--mapped-reads <filename>', String, 'A file with mapped reads for filtering of singletons (unpaired, but mapped reads).') {|m| options[:mapped] = m}
end

optparse.parse!(ARGV)
optparse.parse!(ARGV << '-h') if options.length < 6

# local
unmapped_reads = options[:unmapped]
input_format = options[:format]
phred_quality = options[:phred]
sequencing_type = options[:sequencing]
prefix = options[:prefix]
anchor_length = options[:anchor]
bowtie_version = options[:bowtie_version]
bowtie_index = options[:bowtie_index]
mapped_reads = options[:mapped]

# global
$logfile = File.open("#{prefix}_preparation.log", 'w')

$logfile.puts "Call: ncSplice_prepareUnmapped.rb -u #{unmapped_reads} -f #{input_format} -p #{prefix} -q #{phred_quality} --sequencing-type #{sequencing_type} -a #{anchor_length} -b #{bowtie_version} -x #{bowtie_index} -m #{mapped_reads}"
$logfile.puts ""

# run
##########################################################################################

begin

	# if fasta, prepare reads directly in fasta-format
	if input_format == 'fasta'
			Analysis.fasta2anchors(unmapped_reads, anchor_length, sequencing_type, "#{prefix}_anchors.fasta")
			Analysis.bowtie_map1(bowtie_version, bowtie_index, "#{prefix}_anchors.fasta", input_format, "#{prefix}.bam")	
	
	else
		# if bam, transform unmapped reads into fastq-format before preparing anchors
		if input_format == 'bam'
			Open3.popen3("samtools view #{unmapped_reads}") do |stdin, stdout, stderr, t|
				Analysis.bam2fastq(stdout, "#{prefix}_unmapped.fastq", phred_quality)
			end
			Analysis.prepare_anchorpairs("#{prefix}_unmapped.fastq", anchor_length, sequencing_type, "#{prefix}_anchors.fastq")	
		
		# if fastq, prepare anchors from fastq-format
		else
			Analysis.prepare_anchorpairs(unmapped_reads, anchor_length, sequencing_type, "#{prefix}_anchors.fastq")
		end	

		# map in fastq format
		Analysis.bowtie_map1(bowtie_version, bowtie_index, "#{prefix}_anchors.fastq", input_format, "#{prefix}.bam")	
	end
		
	if sequencing_type == 'pe'
		Open3.popen3("samtools view -bhf 8 #{mapped_reads} > singletons.bam")
	end
	
rescue StandardError => err
	$logfile.puts "#{Time.new}: Error in ncSplice_prepareUnmapped.rb"
	$logfile.puts err.message
	err.backtrace.each {|line| $logfile.puts line}
	exit
end
