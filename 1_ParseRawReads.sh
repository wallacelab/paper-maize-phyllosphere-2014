#! /bin/bash

# Perform minimum entropy decomposition on the Hiseq reads

rawdir=$1
scriptdir=$2
parsedir=$3

join_commands=$parsedir/1_join_commands.txt
echo "echo 'Joining paired ends'" > $join_commands
i=1
for read1 in $rawdir/*_R1.fastq.gz; do
  read2=${read1/_R1.fastq.gz/_R2.fastq.gz}
  outdir=${read1/$rawdir\//}
  outdir=$parsedir/${outdir/.fastq.gz/}

  # if [ ! -e $outdir ]; then continue; fi # Uncomment this if want to continue an interrupted analysis without overwriting previously finished files
  echo "join_paired_ends.py --pe_join_method fastq-join -f $read1 -r $read2 -o $outdir  --perc_max_diff 15 --min_overlap 40 && gzip -f $outdir/*.fastq " >> $join_commands  # Run end-joining then gzip the results
  
  i=$(($i+1))
done
cat $join_commands | parallel --progress
echo -e "\tFinished pair joining!"


# Add sample names to fasta header
commands=$parsedir/1a_sample_commands.txt
echo "echo 'Adding padding and sample names'" > $commands
for readdir in $parsedir/*/; do
  readdir=${readdir%/}
  sample=`echo $readdir | sed -r "s/.+_(LMA._.+_14A.+)_........_........_R./\1/" `  # convoluted way to extract sample name from folder
  echo " zcat $readdir/fastqjoin.join.fastq.gz | fastq_to_fasta | sed \"s/>/>$sample|/\"  > $readdir/joined_with_sample_name.fasta &&  gzip -f $readdir/joined_with_sample_name.fasta " >> $commands
done
cat $commands | parallel --progress
echo -e "\tFinished adding samples!"


# Get read counts
readcount=$parsedir/1b_read_counts.txt
echo -e "sample\treads" > $readcount
for readdir in $parsedir/*/; do
  readdir=${readdir%/}
  sample=`echo $readdir | sed -r "s/.+_(LMA._.+_14A.+)_........_........_R./\1/" `  # convoluted way to extract sample name from folder
  count=$(zcat $readdir/joined_with_sample_name.fasta.gz | wc -l | cut -f1)
  count=$(($count/2))   # 2 lines per fasta entry
  echo -e "$sample\t$count" >> $readcount
done