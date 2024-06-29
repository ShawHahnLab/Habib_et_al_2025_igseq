#!/usr/bin/env bash
cut -f 1 metadata/biosamples.tsv  | tail -n +2 | while read samp; do
	old=$(echo /home/jesse/penn/projects/rmigseq/analysis/merge/*/${samp}.fastq.gz)
	new=analysis/merge/${samp}.fastq.gz
	if [[ -e "$old" && -e "$new" ]]; then
		# skip gzip header and just compare the data
		old_check=$(tail -c +11 $old | md5sum | cut -f 1 -d ' ')
		new_check=$(tail -c +11 $new | md5sum | cut -f 1 -d ' ')
		checks=$(echo -e "$old_check\n$new_check" | uniq)
		echo $samp $checks
	else
		echo $samp ?
	fi
done
