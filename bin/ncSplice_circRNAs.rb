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
require_relative "../lib/circularSplicing.rb"

# options
##########################################################################################

options = {}
optparse = OptionParser.new do |opts|
  opts.banner = "ncSplice version 0.0.0 by Franziska Gruhl (franziska.gruhl@swiss.ch)\n
  Usage: ncSplice_circRNAs.rb -u <unmapped reads> -f <file format: fasta or fastq> -p <prefix> --sequencing-type <se or pe> --library-type <fr-unstranded, fr-firststrand or fd-secondstrand> -l <read length> -b <bowtie-version> -c <chromosomes> [further options]"

  opts.on('-h', '--help', 'Display help screen.') do
    puts opts
    exit
  end
	
	opts.on('-v', '--version', 'Print ncSplice version and dependencies.') do
		puts '# ncSplice'
		puts 'ncSplice v0.1.0'
		puts "\n# Dependencies"
		puts ['ruby >= 2.0.0', 'samtools >= 1.0.0', 'bowtie2 >= 2.1.0 || bowtie >= 0.12.9'].join("\n")
    exit
  end

  options[:anchor] = 20
  options[:bowtie_version] = 'bowtie'
	
 	opts.on('-u', '--unmapped <filename>', String, 'An input file with unmapped reads in fastq format. This file should exist after running ncSplice_prepareUnmapped.rb.') {|u| options[:unmapped] = u}
 	opts.on('-f', '--format <string>', String, 'The input format of unmapped reads file, can be fastq or fasta.') {|f| options[:format] = f}
 	opts.on('--sequencing-type <string>', String, 'The sequencing type. Use \'se\' for single-end and \'pe\' for paired-end sequencing.') {|seq| options[:sequencing] = seq}
	opts.on('--library-type <string>', String, 'The library type. This can be \'fr-unstranded\' for unstranded libraries, and \'fr-firststrand\' or \'fr-secondstrand\' for stranded library types.') {|lib| options[:library] = lib}
  opts.on('-p', '--prefix <string>', String, 'A prefix or label, which will be used for all output files.') {|p| options[:prefix] = p}
  opts.on('-a', '--anchor-length <integer>', Integer, 'The length of the read anchor used in the remapping procedue, the default is 20 bp. Shorter anchors will decrease the mapping precision and longer anchors will cause a reduction in the number of candidates. Make sure that ncSplice_prepareUnampped.rb and ncSplice_circRNAs are supplied with the same anchor length.') {|a| options[:anchor] = a}
  opts.on('-l', '--read-length <integer>', Integer, 'The length of the sequencing read.') {|l| options[:readlength] = l}  
  opts.on('-c', '--chromosomes <directory>', String, 'A directory, which contains all chromosomes to which the unmapped reads were mapped. Sequences need to be in fasta format and in seperated files (one fasta-file/chromosome).') {|c| options[:chromosomes] = c}
  opts.on('-b', '--bowtie-version <String>', String, 'The bowtie version to use, can be either bowtie or bowtie2 (default bowtie).') {|b| options[:bowtie_version] = b}
  opts.on('-s', '--mapped-reads <filename>', String, 'A file with mapped reads for filtering of singletons (unpaired, but mapped reads).') {|s| options[:mapped] = s}
end

optparse.parse!(ARGV)
optparse.parse!(ARGV << '-h') if options.length < 8

# local
unmapped_reads = options[:unmapped]
input_format = options[:format]
sequencing_type = options[:sequencing]
library_type = options[:library]
prefix = options[:prefix]
anchor_length = options[:anchor]
read_length = options[:readlength]
options[:chromosomes][-1] == '/' ? chromosome_files = options[:chromosomes] : chromosome_files = "#{options[:chromosomes]}/"
bowtie_version = options[:bowtie_version]
bowtie_version == 'bowtie' ? bt = '*.ebwt' : bt = '*.bt2'
singletons = options[:mapped]

# global
$logfile = File.open("#{prefix}_mapping.log", 'w')

$logfile.puts "Call: ncSplice_circRNAs.rb -u #{unmapped_reads} -f #{input_format} -p #{prefix} --sequencing-type #{sequencing_type} --library-type #{library_type} -a #{anchor_length} -l #{read_length} -b #{bowtie_version} -c #{chromosome_files} -s #{singletons}"
$logfile.puts ""


# run
##########################################################################################

begin	
 	anchor_pairs = Open3.popen3("samtools view #{prefix}.bam") do |stdin, stdout, stderr, t|
 		CircRNA.process_bam(stdout, chromosome_files)
 	end
 
 	CircRNA.seed_extension(anchor_pairs, anchor_length, read_length, chromosome_files, sequencing_type, library_type, "#{prefix}_candidateReads.txt")	
 	CircRNA.collaps_qnames("#{prefix}_candidateReads.txt", "#{prefix}_candidates.txt")
 	CircRNA.candidates2fa("#{prefix}_candidates.txt", chromosome_files, read_length, "#{prefix}_faIndex.fa")
 	Analysis.bowtie_build(bowtie_version, "#{prefix}_faIndex.fa", "#{prefix}")
 	
 	Open3.popen3("mkdir #{prefix}_index") if !Dir.exists?("#{prefix}_index")
 	Open3.popen3("mv #{prefix}_faIndex.fa #{prefix}_index/; mv #{bt} #{prefix}_index/")
 	
 	Analysis.bowtie_map2(bowtie_version, "#{prefix}_index/#{prefix}", unmapped_reads, input_format, "#{prefix}_remapping.bam")
 
 	Open3.popen3("samtools view #{prefix}_remapping.bam") do |stdin, stdout, stderr, t|
 		CircRNA.remapped_reads(stdout, "#{prefix}_remapped.txt", singletons, read_length, sequencing_type)
 	end
 	
 	CircRNA.final_candidates("#{prefix}_remapped.txt", "#{prefix}_final.txt")
	
rescue StandardError => err
	$logfile.puts "#{Time.new}: Error in ncSplice_circRNAs.rb"
	$logfile.puts err.message
	err.backtrace.each {|line| $logfile.puts line}
	exit
end
