#!/bin/bash

usage()
{
	echo -e "Usage: $0 -d <sample_dir> -p <sample_prefix> -m <modality>"
	exit 1
}

checkDir()
{
	if [ ! -d "$1" ]; then
		echo "Error: $1 directory not found"
		exit 1
	fi
}

checkFile()
{
	if [ ! -f "$1" ]; then
		echo "Error: File $1 not found"
		return 0
	else
		return 1
	fi
}

while getopts "d:p:m:" opts
do
	case "${opts}" in
		d ) sample_dir=$OPTARG ;;
		p ) sample_pre=$OPTARG ;;
		m ) modality=$OPTARG ;;
		? ) usage ;;
	esac
done

checkDir "${sample_dir}"
samples=($(ls "${sample_dir}" | grep "^${sample_pre}"))

for i in "${!samples[@]}"
do
	samples[$i]="${sample_dir}/${samples[$i]}"
	checkDir "${samples[$i]}"
	sample_mod=$(ls "${samples[$i]}" | grep ".*-${modality}-.*")
	samples[$i]="${samples[$i]}/$sample_mod"
	checkDir "${samples[$i]}"
	sample_bam=$(ls "${samples[$i]}" | grep .*bam$)
	samples[$i]="${samples[$i]}/$sample_bam"
	checkFile "${samples[$i]}"
	if [ $? -eq 0 ]; then
		unset 'samples[$i]'
	fi
done

samples=("${samples[@]}")

printf '%s\n' "${samples[@]}"
