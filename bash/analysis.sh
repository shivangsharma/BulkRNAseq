#!/bin/bash

raw_data_path="/media/singhlab/B338-CE9D/Cell_line_RNAseq/data/raw"

if [ -d "$raw_data_directory" ] then
	cd "$raw_data_directory" || exit
	samples=$(find . -type f | sed 's/\.\///' | cut -d'_' -f1-3 | uniq)
else
	echo "Error: Not a directory"
fi

seq_5p=()
seq_3p=()

# /#/ helps prepend a string to all elements of array 
# seq_5p=( "${samples[@]/%/_R1_001.fastq.gz}" )
# seq_3p=( "${samples[@]/%/_R2_001.fastq.gz}" )

for ((i=0; i<${#samples[@]}; i++)); do
	seq_5p[$i]="${samples[$i]_R1_001.fastq.gz}"
	seq_3p[$i]="${samples[$i]_R2_001.fastq.gz}"
done

for ((i=0; i<${#samples[@]}; i++))
do
	cd "$raw_data_path"
	
	fastqc \
	-o /media/singhlab/B338-CE9D/Cell_line_RNAseq/QC/pretrim \
	-dir /media/singhlab/B338-CE9D/Cell_line_RNAseq/tmp \
	-t 40 \
	--nogroup \
	"${seq_5p[i]}"
	
	fastqc \
	-o /media/singhlab/B338-CE9D/Cell_line_RNAseq/QC/pretrim \
	-dir /media/singhlab/B338-CE9D/Cell_line_RNAseq/tmp \
	-t 40 \
	--nogroup 
	"${seq_3p[i]}"
	
	skewer \
	-x AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-m pe \
	-t 32 \
	-o ../skewer-trimmed/23257XR-01-01_S208_L005 \
	-z \
	"${seq_5p[$i]}" "${seq_3p[$i]}"
	
	trimmed_data_path=$(cd ../skewer-trimmed | pwd)
	
	# Need new seq_3p and seq_5p arrays
	
	fastqc \
	-o /media/singhlab/B338-CE9D/Cell_line_RNAseq/QC/posttrim \
	-dir /media/singhlab/B338-CE9D/Cell_line_RNAseq/tmp \
	-t 40 \
	--nogroup \
	"${seq_5p[$i]}"
	
	fastqc \
	-o /media/singhlab/B338-CE9D/Cell_line_RNAseq/QC/posttrim \
	-dir /media/singhlab/B338-CE9D/Cell_line_RNAseq/tmp \
	-t 40 \
	--nogroup 
	"${seq_3p[$i]}"
	
done

axel -n 40 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz -o /home/singhlab/RNA-seq_tools/references
axel -n 20 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz -o /home/singhlab/RNA-seq_tools/references
axel -n 20 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gff3.gz -o /home/singhlab/RNA-seq_tools/references

~/RNA-seq_tools/STAR/STAR-2.7.11a/bin/Linux_x86_64/STAR \
--runThreadN 30 \
--runMode genomeGenerate \
--genomeDir /home/singhlab/RNA-seq_tools/genome_build_STAR \
--genomeFastaFiles /home/singhlab/RNA-seq_tools/references/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /home/singhlab/RNA-seq_tools/references/gencode.v45.primary_assembly.annotation.gtf \
--sjdbOverhang 149

# Pick files that end with fastq.gz
skewer_data_dir="/media/singhlab/B338-CE9D/Cell_line_RNAseq/data/skewer-trimmed/"
if [ -d "$skewer_data_dir" ]; then
	cd "$skewer_data_dir"
	samples=$( find . -type f -name "*.fastq.gz" | sed 's/\.\///' | cut -d'-' -f1-4 | uniq )
	readarray -t samples <<< "$samples"
else
	echo "Error: Not a directory"
fi

# seq_5p=("${samples[@]/%/-pair1.fastq.gz}")
# seq_3p=("${samples[@]/%/-pair2.fastq.gz}")
seq_5p=()
seq_3p=()
for ((i=0; i<${#samples[@]}; i++)); do
	seq_5p[$i]="$skewer_data_dir${samples[$i]}-pair1.fastq.gz"
	seq_3p[$i]="$skewer_data_dir${samples[$i]}-pair2.fastq.gz"
done

join () {
  local IFS="$1"
  shift
  echo "$*"
}

# file_list_5p=$(join "," "${seq_5p[@]}")
# file_list_3p=$(join "," "${seq_3p[@]}")

for ((i=0; i<${#samples[@]}; i++)); do
	~/RNA-seq_tools/STAR/STAR-2.7.11a/bin/Linux_x86_64/STAR \
	--runThreadN 30 \
	--genomeDir /home/singhlab/RNA-seq_tools/genome_build_STAR \
	--readFilesIn "${seq_5p[$i]}" "${seq_3p[$i]}" \
	--readFilesCommand zcat \
	--outFileNamePrefix "/media/singhlab/B338-CE9D/Cell_line_RNAseq/data/STAR_alignment/${samples[$i]}" \
	--outTmpDir /home/singhlab/RNA-seq_tools/tmp/ \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--chimSegmentMin 18 \
	--quantMode TranscriptomeSAM \
	--sjdbOverhang 149
	--twopass1readsN -1
done

~/RNA-seq_tools/RSEM/rsem-prepare-reference \
--gtf ~/RNA-seq_tools/references/gencode.v45.primary_assembly.annotation.gtf \
~/RNA-seq_tools/references/GRCh38.primary_assembly.genome.fa \
~/RNA-seq_tools/transcriptome_reference_RSEM/RSEM_ref

for ((i=0; i<${#samples[@]}; i++)); do
	~/RNA-seq_tools/RSEM/rsem-calculate-expression \
	--bam \
	--no-bam-output \
	--num-threads 30 \
	--paired-end \
	--forward-prob 0 \
	"/media/singhlab/B338-CE9D/Cell_line_RNAseq/data/STAR_alignment/${samples[$i]}Aligned.toTranscriptome.out.bam" \
	~/RNA-seq_tools/transcriptome_reference_RSEM/RSEM_ref \
	"/media/singhlab/B338-CE9D/Cell_line_RNAseq/data/RSEM_quantification/${samples[$i]}"
done
