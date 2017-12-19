#! /bin/bash

# Check HiSeq runs

outdir=$1
rawdir=$2
scriptdir=$3
plate_key=$4
md5_sums=$5

# Check MD5 sums 
echo "Checking MD5 sums of files"
md5file=$outdir/0a_round2_hiseq_md5s.downloaded.txt
echo "file	md5sum" > $md5file
for file in $rawdir/*.gz; do
  name=${file/*\//}	# Strip leading directory
  md5=$(md5sum $file | cut -f1 -d" ")
  echo "$name	$md5" >> $md5file
done
python3 $scriptdir/0b_CompareMd5s.py -s $md5_sums -d $md5file -o $outdir/0b_round2_hiseq_md5s.mismatches.txt


# Get raw read counts and confirm blanks
echo "Tabulating raw read counts and confirming blanks"
readcount=$outdir/0c_raw_read_counts.txt
echo "sample	reads" > $readcount
for infile in $rawdir/*R1.fastq.gz; do	# Only take forward read
  reads=$(zcat $infile | wc -l)
  reads=$(($reads / 4))
  sample=$(echo $infile | sed -r "s|.+_HJC7WADXX_||" | sed -r "s|_[ACGT]+_[ACGT]+_R[12].fastq.gz||")
  echo $sample
  echo "$sample	$reads">>$readcount
done
Rscript $scriptdir/0d_DoReadCountQC.r $readcount $outdir/0d_read_count_qc.png
python3 $scriptdir/0d_PlotReadCountsByPlates.py -i $readcount -k $plate_key -o $outdir/0d_reads_by_plates.png

