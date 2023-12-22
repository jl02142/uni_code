# J ludwig 2017

#PBS -S /bin/bash
#PBS -q batch
#PBS -N bwa_mem_generate
#PBS -l nodes=1:ppn=1:dbnode
#PBS -l walltime=24:00:00
#PBS -l mem=10gb

#PBS -m ae
for sample in *_R1.fastq.gz
do
# fix dir path below *NAME*
if [ ! -e /lustre1/NAME/Candida/new_seq/temp_jobid ]; then echo $sample >> /lustre1/NAME/Candida/new_seq/temp_jobid; else echo $sample >> /lustre1/NAME/Candida/new_seq/temp_jobid; fi
echo "#PBS -S /bin/bash" > $sample.del.sh
echo "#PBS -q batch" >> $sample.del.sh
echo "#PBS -N $sample" >> $sample.del.sh
echo "#PBS -l nodes=1:ppn=1:dbnode" >> $sample.del.sh
echo "#PBS -l walltime=24:00:00" >> $sample.del.sh
echo "#PBS -l mem=10gb" >> $sample.del.sh
echo "#PBS -m ae" >> $sample.del.sh
echo "cd /lustre1/NAME/Candida/new_seq" >> $sample.del.sh
echo "module load bwa/0.7.10" >> $sample.del.sh
echo "module load samtools/1.2" >> $sample.del.sh
echo "sample=$sample" >> $sample.del.sh
echo "input1=$sample" >> $sample.del.sh
echo 'input2=`basename $sample _R1.fastq.gz`_R2.fastq.gz' >> $sample.del.sh
echo 'output=`basename $sample _R1.fastq.gz`' >> $sample.del.sh
echo 'bwa mem -c 20000 C_albicans_SC5314_A.fa $input1 $input2   | samtools view -b - | samtools sort - /lustre1/NAME/Candida/new_seq/analysis/BWA_c_20000/$output' >> $sample.del.sh
list=${sample%%_L00*_R*.fastq.gz}
list2=${list#Run*_}
if [ ! -e /lustre1/NAME/Candida/new_seq/list_id ]; then grep -q -F "$list2" list_id || echo "$list2" >> list_id; else grep -q -F "$list2" list_id || echo "$list2" >> list_id; fi
if [ ! -e /lustre1/NAME/Candida/new_seq/temp_jobid2 ]; then echo $list2 >> /lustre1/NAME/Candida/new_seq/temp_jobid2; elif grep -q "$list2" temp_jobid2 ; then :; else echo $list2 >> /lustre1/NAME/Candida/new_seq/temp_jobid2; 
	echo "#PBS -S /bin/bash" > $list2.del2.sh
	echo "#PBS -q batch" >> $list2.del2.sh
	echo "#PBS -N $sample 2" >> $list2.del2.sh
	echo "#PBS -l nodes=1:ppn=1:dbnode" >> $list2.del2.sh
	echo "#PBS -l walltime=24:00:00" >> $list2.del2.sh
	echo "#PBS -l mem=10gb" >> $list2.del2.sh
	echo "#PBS -m ae" >> $list2.del2.sh
	echo "cd /lustre1/NAME/Candida/new_seq/analysis/" >> $list2.del2.sh
	echo "module load samtools/1.2" >> $list2.del2.sh
	i=1
	for thing2 in *$list2*.gz
	do
	if [[ $thing2 =~ R1.fastq.gz ]]; then thing3=${thing2%%_R1.fastq.gz}
	echo "input"$i"="$thing3.bam >> $list2.del2.sh
	input[i]="$"input$i
	((i++))
	fi
	done
	echo "output="$list2 >> $list2.del2.sh
	echo 'samtools merge /lustre1/NAME/Candida/new_seq/analysis/BWA_c_20000/$output.bam' ${input[*]} >> $list2.del2.sh
	fi
done
while read input; do qcheck=$(qsub $input.del.sh);temp1=${input%_[A-Z]*-[A-Z]*_L00*_R*.fastq.gz};temp2=${temp1#Run*_};qcheck2=${qcheck%%.scm}; echo "$temp2 $qcheck2" >> check3; eval "$temp2+=$qcheck2:";  done < temp_jobid
while read input2; do input3=${input2%_[A-Z]*-[A-Z]*};echo "$input2 $input3" >> check4; qsub -W depend=afterok:${!input3} $input2.del2.sh; done < list_id
rm temp_jobid
rm list_id
rm temp_jobid2
