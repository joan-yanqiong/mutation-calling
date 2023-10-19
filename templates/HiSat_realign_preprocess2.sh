#!/bin/sh

echo "Loading modules..."
module load samtools/1.14
module load bamtools/2.4.2
module load plinkseq/0.10

echo "Creating output directory..."
mkdir -p !{sample_id}
cd !{sample_id}

echo "Setting up variables..."
LIB="../!{jar_files}"
MAF="../!{maf}"
BAM="../!{recal_bam}"
REFo="../!{ref_path}"
NAME=!{sample_id}

xbase=${BAM##*/}
STUB=${xbase%.*}


echo "LIB: $LIB"
echo "MAF: $MAF"
echo "BAM: $BAM"
echo "REFo: $REFo"
echo "NAME: $NAME"
echo "STUB: $STUB"

#WORK_DIR=$(echo pwd) #working directory

#check to see if maflite or maf file:
#maf lite will have columns: contig, start_position, end_position
#maf will have columns: Chromosome Start_position End_position

##########################################################
#Generate Intervals
##########################################################

echo "generate intervals from maf or maflite"

titlemaf=`cat $MAF | head -n100 | grep -P "position\s+" | sed 's/\t/\n/g' | grep -P ^Chromosome`
title=`cat $MAF | head -n100 | grep -P "chr\s+" | sed 's/\t/\n/g' | grep -P ^chr`
if [ ! -z "$titlemaf" ];
then
echo "this is a maf file"
grep -v "#" $MAF | awk -F '\t' '{print $5, $6, $7}' | sed 's/ /\t/g' | sed -e 's/\s\+/:/'| sed -e 's/\s\+/-/' | sed '/^#/ d' | sed -e '1d' | sed 's/^M:/MT:/' > snp_mutations.intervals
# grep -v "#" $MAF | gcol Chromosome Start_position End_position | sed -e 's/\s\+/:/'| sed -e 's/\s\+/-/' | sed '/^#/ d' | sed -e '1d' | sed -e 's/^M:/MT:/' > snp_mutations.intervals
if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi
elif [ ! -z "$title" ];
then
echo "this is a maflite"
grep -v "#" $MAF | awk -F '\t' '{print $5, $6, $7}' | sed 's/ /\t/g' | sed -e 's/\s\+/:/'| sed -e 's/\s\+/-/' | sed '/^#/ d' | sed -e '1d' | sed 's/^M:/MT:/' > snp_mutations.intervals
# grep -v "#" $MAF | gcol chr start end | sed -e 's/\s\+/:/'| sed -e 's/\s\+/-/' | sed '/^#/ d' | sed -e '1d' | sed 's/^M:/MT:/' > snp_mutations.intervals
if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi

fi

#print out intervals
sed 's/:/'"$(printf '\011')"'/g' snp_mutations.intervals | sed 's/-/'"$(printf '\011')"'/g' > snp_mutations.intervals.bed

#######################################################
#Extracting Paired Reads from BAM
#######################################################

echo "extracting paired reads from bam"

samtools view -L snp_mutations.intervals.bed $BAM | cut -f1 > IDs_all.txt
if [ "$?" -ne 0 ]; then echo "command extracting read IDs failed"; exit 1; fi

java -Xmx7g -jar $LIB/FilterSamReads.jar I=$BAM O=tmp_bam.bam READ_LIST_FILE=IDs_all.txt FILTER=includeReadList WRITE_READS_FILES=false VALIDATION_STRINGENCY=LENIENT
if [ "$?" -ne 0 ]; then echo "command generating BAM based on read list failed"; exit 1; fi


echo "convert bam to fastq"
samtools view -H tmp_bam.bam | sed '$d' - > tmp_header_T.sam
# if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi

samtools view tmp_bam.bam | awk '$2 < 2040 { print }' > tmp0_T.sam
# if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi

cat tmp_header_T.sam tmp0_T.sam > tmp_filteredbamT.sam

# Added || true, to suppress/ignore the error.
java -Xmx7g -jar $LIB/SamToFastq.jar I=tmp_filteredbamT.sam F=${NAME}_tmp_sequence_1.fastq F2=${NAME}_tmp_sequence_2.fastq VALIDATION_STRINGENCY=LENIENT || true
# if [ "$?" -ne 0 ]; then echo "command failed"; exit 1; fi
