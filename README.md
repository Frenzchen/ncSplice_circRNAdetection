# ncSplice_circRNAdetection
This repository holds the part from ncSplice that was used in my PhD thesis to detect circRNAs from RNA sequencing data

## Filtering of unmapped reads from Illumina RNA sequencing data.
The file unmapped.bam as provided by TopHat2 was used as input for ncSplice. This part of the pipeline will filter the 
unmapped reads for good-quality reads, prepare a fastq-file with the anchor pairs and map them on the provided index, if
the option -m is provided, it will scan the accepted_hits.bam for singletons. Latter are used to filter backsplice junctions
(one part of the read is on the backsplice junction and the second part should be inside the circRNA).

Command as used in the thesis:

`ncSplice_prepareUnmapped.rb -u unmapped.bam -f bam -p myFile --sequencing-type pe -a 20 -b bowtie2 -x bt2-index -m accepted_hits.bam`

## circRNA detection
This part will perform the seed extension of mapped anchors, create a new index from the backslice junctions to remap all unmapped reads
and summarize results in a final output file

Command as used in the thesis:

`ncSplice_circRNAs.rb -u myFile_unmapped.fastq -f fastq -p myFile --sequencing-type pe --library-type fr-firststrand -a 20 -l 101 -b bowtie2 -c chromosomes/ -s singletons.bam`

Note: In the next step, candidates will be compared to RNase R treated samples to filter the putative circRNA candidates. If you don't
have RNase R treated samples, it is not recommended to use this set of scripts, because the pipeline was optimized for comparison.
