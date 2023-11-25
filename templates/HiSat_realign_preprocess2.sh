#!/bin/sh

echo "$(date)   Loading modules..."
module load samtools/1.14
module load bamtools/2.4.2
module load plinkseq/0.10

echo "$(date)  Creating output directory..."
mkdir -p !{sample_id}
cd !{sample_id}

echo "$(date)  Setting up variables..."
LIB="../!{jar_files}"
MAF="../!{maf}"
BAM="../!{recal_bam}"
REFo="../!{ref_path}"
NAME="!{sample_id}"

echo "$(date)   Determine prefix for output files..."
if [[ "!{sample_type}" == "tumor" ]]
then
    prefix="!{sample_id}"
elif [[ "!{sample_type}" == "normal" ]]
then
    prefix="!{pair_id}"
fi

xbase=${BAM##*/}
STUB=${xbase%.*}


echo "LIB: $LIB"
echo "MAF: $MAF"
echo "BAM: $BAM"
echo "REFo: $REFo"
echo "NAME: $NAME"
echo "STUB: $STUB"

#check to see if maflite or maf file:
#maf lite will have columns: contig, start_position, end_position
#maf will have columns: Chromosome Start_position End_position

################################
# ---- Generate Intervals ---- #
################################

echo "$(date)   Generate intervals from maf or maflite and save..."

titlemaf=`cat $MAF | head -n100 | grep -P "position\s+" | sed 's/\t/\n/g' | grep -P ^Chromosome`
title=`cat $MAF | head -n100 | grep -P "chr\s+" | sed 's/\t/\n/g' | grep -P ^chr`

if [ ! -z "$titlemaf" ];
then
echo "$(date)   Type of file: maf"
grep -v "#" $MAF | awk -F '\t' '{print $5, $6, $7}' | sed 's/ /\t/g' | sed -e 's/\s\+/:/'| sed -e 's/\s\+/-/' | sed '/^#/ d' | sed -e '1d' | sed 's/^M:/MT:/' > snp_mutations.intervals

# Dealing with maflite files
if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi
elif [ ! -z "$title" ];
then
echo "$(date)   Type of file: maflite"
grep -v "#" $MAF | awk -F '\t' '{print $5, $6, $7}' | sed 's/ /\t/g' | sed -e 's/\s\+/:/'| sed -e 's/\s\+/-/' | sed '/^#/ d' | sed -e '1d' | sed 's/^M:/MT:/' > snp_mutations.intervals

if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi

fi

echo "$(date)   Save intervals to bed file..."
sed 's/:/'"$(printf '\011')"'/g' snp_mutations.intervals | sed 's/-/'"$(printf '\011')"'/g' > snp_mutations.intervals.bed

##############################################
# ---- Extracting Paired Reads from BAM ---- #
##############################################

echo "$(date)   Extracting paired reads from BAM..."

samtools view -L snp_mutations.intervals.bed $BAM | cut -f1 > IDs_all.txt
if [ "$?" -ne 0 ]; then echo "Extracting read IDs failed"; exit 1; fi

java -Xmx7g -jar $LIB/FilterSamReads.jar I=$BAM O=tmp_bam.bam READ_LIST_FILE=IDs_all.txt FILTER=includeReadList WRITE_READS_FILES=false VALIDATION_STRINGENCY=LENIENT
if [ "$?" -ne 0 ]; then echo "Generating BAM based on read list failed"; exit 1; fi

echo "$(date)   Convert BAM to FASTQ..."
samtools view -H tmp_bam.bam | sed '$d' - > tmp_header_T.sam

samtools view tmp_bam.bam | awk '$2 < 2040 { print }' > tmp0_T.sam

cat tmp_header_T.sam tmp0_T.sam > tmp_filteredbamT.sam

# IMPORTANT: Added || true, to suppress/ignore the error.
java -Xmx7g -jar $LIB/SamToFastq.jar I=tmp_filteredbamT.sam F=${prefix}_tmp_sequence_1.fastq F2=${prefix}_tmp_sequence_2.fastq VALIDATION_STRINGENCY=LENIENT || true

touch "ok.txt"

echo "$(date)   Completed"