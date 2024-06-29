#!/usr/bin/env bash

# All match up except 44715 due to a tiny extra sample I included in this
# version (and the output matches, so effectively all the inputs are the same)

ref=(kimdb sonarramesh sonarramesh)
there=(mu kappa lambda)
here=(IGH IGK IGL)
grep IgM metadata/biosamples.tsv | cut -f 16 | sort -u | while read subject; do
	for idx in $(seq 0 2); do
		locus=${here[idx]}
		old=$(echo /home/jesse/penn/projects/rmigseq/analysis/igdiscover/${ref[idx]}/${there[idx]}/${subject}/reads.fastq.gz)
		new=analysis/igdiscover/${ref[idx]}/${locus}/${subject}/reads.fastq.gz
		if [[ -e "$old" && -e "$new" ]]; then
			# skip gzip header and just compare the data
			# (NOPE this won't work for cases where we've just literally
			# cat'd gzip files together, like is the case for some
			# of our IgDiscover input files.)
			#old_check=$(tail -c +11 $old | md5sum | cut -f 1 -d ' ')
			#new_check=$(tail -c +11 $new | md5sum | cut -f 1 -d ' ')
			old_check=$(zcat $old | md5sum | cut -f 1 -d ' ')
			new_check=$(zcat $new | md5sum | cut -f 1 -d ' ')
			checks=$(echo -e "$old_check\n$new_check" | uniq)
			echo $subject $locus $checks
		else
			echo $subject $locus ?
		fi
	done
done
