#########################################################################################################################
#	This shell script takes input as folder name of Normal and Tumor followed by list of files of Normal and Tumor.	#
#	e.g WGS_pipeline_TumorVsNormal_hg38_v0.1.sh Blood Tumor list_Blood list_Tumor					#
#	Last Modified on 29/12/2022											#
#	Last modified: Adding vcf reformat.										#
#															#
#															#
#															#
#	Tools Used in this pipeline											#
#	1.  Fastqc													#
#	2.  Trimmomatic													#
#	3.  BWA														#
#	4.  Markduplicate												#
#	5.  BQSR													#
#	6.  HaplotypeCaller (Germline)											#
#	7.  Strelka Germline (Germline)											#
#	8.  Platypus (Germline)												#
#	9.  Mutect2 (Somatic)												#
#	10. Varscan (Somatic)												#
#	11. Lofreq (Somatic)												#
#	12. CNVpytor (For copy number variation)									#
#	13. Genefuse (For Fusion gene analysis)										#
#	14. MSI sensor (For Microsatellite Instability Analysis)							#
#	15. Delly and Pindel (For Structural variant analysis)								#
#	16. Annovar (Variant Annotation)										#
#########################################################################################################################

# Taking user inputs and setting resource paths
input1="$1"	# Taking input of normal folder name as argument
input2="$2"	# Taking input of tumor folder name as argument
input3="$3"	# Taking input of list file of normal sample
input4="$4"	# Taking input of list file of tumor sample
workdir=$(pwd)	# storing path of Working directory
raw_normal="$workdir/$input1"	# Assigning path of normal sample
raw_tumor="$workdir/$input2"	# Assigning path of tumor sample 
out="$workdir/Output_hg38_${input1}_Vs_${input2}"	# Assigning path of output folder
db="/mnt/database"		# Assigning database folder path
hg38="$db/GRCH38.P14/hg38.fa"	# Assigning path of human genome reference file


# Removing pre-exiting files if any

rm $out/WGS_analysis.log
rm $out/parallel_fastqc_Normal
rm $out/parallel_trimmomatic_Normal
rm $out/parallel_after_trimmomatic_Normal
rm $out/parallel_fastqc_Tumor
rm $out/parallel_trimmomatic_Tumor
rm $out/parallel_after_trimmomatic_Tumor
rm $out/parallel_BWA
rm $out/parallel_markdup
rm $out/parallel_BaseRecalibrator
rm $out/parallel_ApplyBQSR
rm $out/parallel_mpileup
rm $out/parallel_samtools_indexing
rm $out/parallel_variantCaller
rm $out/parallel_SV

echo "First $input1 Second $input2 Third $input3 Forth $input4" # Print inputs as a test

# Creating folders for the intermediate results
mkdir -p $out/$input1/alignments_stats
mkdir -p $out/$input2/alignments_stats
mkdir -p $out/$input1/bqsr
mkdir -p $out/$input2/bqsr
mkdir -p $out/$input1/markdup
mkdir -p $out/$input2/markdup
mkdir -p $out/$input1/FastQC_output
mkdir -p $out/$input2/FastQC_output
mkdir -p $out/$input1/trimmomatic_output
mkdir -p $out/$input2/trimmomatic_output
mkdir -p $out/$input1/FastQC_output_after-trimmomatic
mkdir -p $out/$input2/FastQC_output_after-trimmomatic
mkdir -p $out/$input1/multiqc
mkdir -p $out/$input2/multiqc
mkdir -p $out/$input1/seqkit
mkdir -p $out/$input2/seqkit


eval "$(conda shell.bash hook)"		# Setting bash for conda environment
conda activate wgs_gatk4		# Activating conda environment

# Fastqc and Trimmomatics run for Normal
while IFS= read -r line			# While loop to add scripts for each lane in "parallel_fastqc_Normal" 
	do 				# "parallel_trimmomatic_Normal" and "parallel_after_trimmomatic_Normal"
		echo $line
		sm=$(echo $raw_normal/$line | xargs -n 1 basename | cut -f 1-3 -d"_"); 	# Extracting file name from file 
		nam=$(echo $raw_normal/$line | xargs -n 1 basename | cut -f 1-2 -d"_");	# list
		echo $sm
		echo $nam

# QC checking

		echo "fastqc $raw_normal/$line -t 30 -o $out/$input1/FastQC_output/" >> $out/parallel_fastqc_Normal

# Adaptor Trimming and cleaning

		echo "trimmomatic PE -threads 30  -phred33 $raw_normal/${sm}_R1* $raw_normal/${sm}_R2* $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R1-trimmed_UP.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_UP.fastq ILLUMINACLIP:$db/Trimmomatic_adaptors/adaptors:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50" >> $out/parallel_trimmomatic_Normal

# QC-rechecking after trimming

		echo "fastqc -t 30 $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq -o $out/$input1/FastQC_output_after-trimmomatic/" >> $out/parallel_after_trimmomatic_Normal

	done < "$input3"

# Running script in Parallel

echo "Fastqc for $input1 Started" >> $out/WGS_analysis.log # Writing in log file
date >> $out/WGS_analysis.log	# Writing in log file

parallel -j 2 < $out/parallel_fastqc_Normal	# Running script in parallel

echo "Fastqc for $input1 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file

echo "Trimmomatic for $input1 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file

parallel -j 2 < $out/parallel_trimmomatic_Normal	# Running script in parallel

echo "Trimmomatic for $input1 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file
echo "##############################" >> $out/WGS_analysis.log		# Writing in log file

echo "Fastqc after trimmomatic for $input1 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file

parallel -j 2 < $out/parallel_after_trimmomatic_Normal		# Running script in parallel

echo "Fastqc after trimmomatic for $input1 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file
echo "##############################" >> $out/WGS_analysis.log			# Writing in log file

# Seqkit for Normal

echo "Seqkit before and after trim for $input1 Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log

seqkit stats -j 60 $raw_normal/* > $out/$input1/seqkit/${nam}.txt
seqkit stats -j 60 $out/$input1/trimmomatic_output/* > $out/$input1/seqkit/${nam}_after_trim.txt

echo "Seqkit before and after trim for $input1 Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log
echo "##############################" >> $out/WGS_analysis.log

# Fastqc and Trimmomatic run for Tumor

while IFS= read -r line
do
	echo $line
	sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");
	nam=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-2 -d"_");
	echo $sm
	echo $nam

# QC checking
		echo "fastqc $raw_tumor/$line -t 30 -o $out/$input2/FastQC_output/" >> $out/parallel_fastqc_Tumor
# Trimming
	
		echo "trimmomatic PE -threads 30 -phred33 $raw_tumor/${sm}_R1* $raw_tumor/${sm}_R2* $out/$input2/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input2/trimmomatic_output/${sm}_R1-trimmed_UP.fastq $out/$input2/trimmomatic_output/${sm}_R2-trimmed_P.fastq $out/$input2/trimmomatic_output/${sm}_R2-trimmed_UP.fastq ILLUMINACLIP:$db/Trimmomatic_adaptors/adaptors:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50" >> $out/parallel_trimmomatic_Tumor

# QC-rechecking-after-trimming

		echo "fastqc -t 30 $out/$input2/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input2/trimmomatic_output/${sm}_R2-trimmed_P.fastq -o $out/$input2/FastQC_output_after-trimmomatic/" >> $out/parallel_after_trimmomatic_Tumor

done < "$input4"

echo "Fastqc for $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file

parallel -j 2 < $out/parallel_fastqc_Tumor		# Running script in parallel

echo "Fastqc for $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file

echo "Trimmomatic for $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file

parallel -j 2 < $out/parallel_trimmomatic_Tumor		# Running script in parallel

echo "Trimmomatic for $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file
echo "##############################" >> $out/WGS_analysis.log		# Writing in log file

echo "Fastqc after trimmomatic for $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file

parallel -j 2 < $out/parallel_after_trimmomatic_Tumor	# Running script in parallel

echo "Fastqc after trimmomatic for $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file
echo "##############################" >> $out/WGS_analysis.log			# Writing in log file


# Seqkit for Tumor sample

echo "Seqkit before and after trim for $input2 Started" >> $out/WGS_analysis.log    # Writing in log file
date >> $out/WGS_analysis.log

seqkit stats -j 60 $raw_tumor/* > $out/$input2/seqkit/${nam}.txt
seqkit stats -j 60 $out/$input2/trimmomatic_output/* > $out/$input2/seqkit/${nam}_after_trim.txt

echo "Seqkit before and after trim for $input2 Started" >> $out/WGS_analysis.log    # Writing in log file
date >> $out/WGS_analysis.log
echo "##############################" >> $out/WGS_analysis.log

# MultiQC before and after trimming

echo "Multiqc for $input1 and $input2 Started" >> $out/WGS_analysis.log    # Writing in log file
date >> $out/WGS_analysis.log

multiqc $out/$input1/FastQC_output/* -n ${input1}_beforetrim_multiqc_report -o $out/$input1/multiqc
multiqc $out/$input2/FastQC_output/* -n ${input2}_beforetrim_multiqc_report -o $out/$input2/multiqc
multiqc $out/$input1/FastQC_output_after-trimmomatic/* -n ${input1}_aftertrim_multiqc_report -o $out/$input1/multiqc
multiqc $out/$input2/FastQC_output_after-trimmomatic/* -n ${input2}_aftertrim_multiqc_report -o $out/$input2/multiqc

echo "Multiqc for $input1 and $input2 Completed" >> $out/WGS_analysis.log    # Writing in log file
date >> $out/WGS_analysis.log
echo "##############################" >> $out/WGS_analysis.log


############# BWA for Normal sample

while IFS= read -r line
do
		
	echo $line
	header1=$(zcat $raw_normal/$line | head -n 1);	# Extracting header from fastq file. Will not work if fastq file
							# is not gunzip file. Change command to cat
	id1=$(echo $header1 | cut -f 3-4 -d":" | sed 's/@//');	# Extracting run ID from fastq file
	sm1=$(echo $raw_normal/$line | xargs -n 1 basename | cut -f 1-3 -d"_");	# Extracting normal file name upto 3 place
										# delimited by '_'
	nam1=$(echo $raw_normal/$line | xargs -n 1 basename | cut -f 1-2 -d"_");  # Extracting normal file name upto
										# 2 places delimited by '_'
	echo "header contant" $header # printing header
	echo $id1
	echo $sm1
	echo $nam1
	echo '@RG\tID:$id\tSM:$sm\tPL:ILLUMINA'
	echo "Read Group @RG\tID:$id1\tSM:$sm1\tPL:ILLUMINA"
	echo "BWA aligner running for" $line
	if [ -z "$(ls -A  $out/$input1/trimmomatic_output/)" ]   # Checking folder if empty
	then
		echo "bwa mem -t 30 -M -R \"@RG\\tID:$id1\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:$sm1\\tPI:200\" $hg38 $raw_normal/${sm1}_R1.fastq $raw_normal/${sm1}_R2.fastq 2> $out/$input1/${sm1}.align.stderr | samtools sort -@ 30 -o $out/$input1/${sm1}.sorted.bam" >> $out/parallel_BWA	# Writing file for parallel run	when trimmomatic_output is empty
	else
		echo "bwa mem -t 30 -M -R \"@RG\\tID:$id1\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:$sm1\\tPI:200\" $hg38 $out/$input1/trimmomatic_output/${sm1}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm1}_R2-trimmed_P.fastq 2> $out/$input1/${sm1}.align.stderr | samtools sort -@ 30 -o $out/$input1/${sm1}.sorted.bam" >> $out/parallel_BWA	# Writing file for 
												# parallel run when
												# trimmmomatic_output is 
												# not empty
	fi
	done < "$input3"	# Passing argument for list normal

# BWA for Tumor sample

while IFS= read -r line					# Same as the above for BWA Normal
do							# for tumor sample
	echo $line
	header2=$(zcat $raw_tumor/$line | head -n 1);
	id2=$(echo $header2 | cut -f 3-4 -d":" | sed 's/@//');
	sm2=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");
	nam2=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-2 -d"_");

	echo "header contant"$header2
	echo $id2
	echo $sm2
	echo $nam2
	echo '@RG\tID:$id\tSM:$sm2\tPL:ILLUMINA'
	echo "Read Group @RG\tID:$id2\tSM:$sm2\tPL:ILLUMINA"
	echo "BWA aligner running for" $line

	if [ -z "$(ls -A  $out/$input2/trimmomatic_output/)" ]
	then
		echo "bwa mem -t 30 -M -R \"@RG\\tID:$id2\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:$sm2\\tPI:200\" $hg38 $raw_tumor/${sm2}_R1.fastq $raw_tumor/${sm2}_R2.fastq 2> $out/$input2/${sm2}.align.stderr | samtools sort -@ 30 -o $out/$input2/${sm2}.sorted.bam" >> $out/parallel_BWA
	else

		echo "bwa mem -t 30 -M -R \"@RG\\tID:$id2\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:$sm2\\tPI:200\" $hg38 $out/$input2/trimmomatic_output/${sm2}_R1-trimmed_P.fastq $out/$input2/trimmomatic_output/${sm2}_R2-trimmed_P.fastq 2> $out/$input2/${sm2}.align.stderr | samtools sort -@ 30 -o $out/$input2/${sm2}.sorted.bam" >> $out/parallel_BWA
	fi

	done < "$input4"

echo "BWA for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file

parallel -j 2 < $out/parallel_BWA	# Running script in parallel

echo "BWA for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file
echo "##############################" >> $out/WGS_analysis.log		# Writing in log file


###MarkDuplicates For Normal

echo "java -XX:+UseParallelGC -XX:ParallelGCThreads=30 -Xms50G -Xmx50G -jar /home/genomics/software/picard-tools-1.119/MarkDuplicates.jar I=$out/$input1/${nam1}_L001.sorted.bam I=$out/$input1/${nam1}_L002.sorted.bam I=$out/$input1/${nam1}_L003.sorted.bam I=$out/$input1/${nam1}_L004.sorted.bam O=$out/$input1/markdup/${nam1}_dedup.bam M=$out/$input1/markdup/${nam1}_metrics.txt 2> $out/$input1/markdup/${nam1}.stderr" >> $out/parallel_markdup		# Writing file for parallel run

###MarkDuplicates For Tumor

echo "java -XX:+UseParallelGC -XX:ParallelGCThreads=30 -Xms50G -Xmx50G -jar /home/genomics/software/picard-tools-1.119/MarkDuplicates.jar I=$out/$input2/${nam2}_L001.sorted.bam I=$out/$input2/${nam2}_L002.sorted.bam I=$out/$input2/${nam2}_L003.sorted.bam I=$out/$input2/${nam2}_L004.sorted.bam O=$out/$input2/markdup/${nam2}_dedup.bam M=$out/$input2/markdup/${nam2}_metrics.txt 2> $out/$input2/markdup/${nam2}.stderr" >> $out/parallel_markdup		# Writing file for parallel run

echo "Markduplicate for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file

parallel -j 2 < $out/parallel_markdup	# Running script in parallel

echo "Markduplicate for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file
echo "##############################" >> $out/WGS_analysis.log			# Writing in log file

echo "Sambamba Alignment Statistics for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log									# Writing in log file

#Alignment Stats For Normal
sambamba flagstat -t 60 $out/$input1/markdup/${nam1}_dedup.bam > $out/$input1/alignments_stats/${nam1}.txt	# Getting 
													# alignment 
													#statistics

samtools index -@ 60 $out/$input1/markdup/${nam1}_dedup.bam	# Generating index from markduplicate bam for normal 
								# sample
#Alignment Stats For Tumor

sambamba flagstat -t 60 $out/$input2/markdup/${nam2}_dedup.bam > $out/$input2/alignments_stats/${nam2}.txt # Generating index 
												# from markduplicate 
												# bam for tumor

samtools index -@ 60 $out/$input2/markdup/${nam2}_dedup.bam 	# Generating index from markduplicate bam for tumor 
								# sample

# Calculating Coverage for Normal Sample for each chromosome

samtools depth -r chr1:0-248956422 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr1_depth_coverage_samtools.txt
samtools depth -r chr2:0-242193529 --reference $hg38  $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr2_depth_coverage_samtools.txt
samtools depth -r chr3:0-198295559 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr3_depth_coverage_samtools.txt
samtools depth -r chr4:0-190214555 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr4_depth_coverage_samtools.txt
samtools depth -r chr5:0-181538259 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr5_depth_coverage_samtools.txt
samtools depth -r chr6:0-170805979 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr6_depth_coverage_samtools.txt
samtools depth -r chr7:0-159345973 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr7_depth_coverage_samtools.txt
samtools depth -r chr8:0-145138636 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr8_depth_coverage_samtools.txt
samtools depth -r chr9:0-138394717 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr9_depth_coverage_samtools.txt
samtools depth -r chr10:0-133797422 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr10_depth_coverage_samtools.txt
samtools depth -r chr11:0-135086622 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr11_depth_coverage_samtools.txt
samtools depth -r chr12:0-133275309 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr12_depth_coverage_samtools.txt
samtools depth -r chr13:0-114364328 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr13_depth_coverage_samtools.txt
samtools depth -r chr14:0-107043718 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr14_depth_coverage_samtools.txt
samtools depth -r chr15:0-101991189 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr15_depth_coverage_samtools.txt
samtools depth -r chr16:0-90338345 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr16_depth_coverage_samtools.txt
samtools depth -r chr17:0-83257441 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr17_depth_coverage_samtools.txt
samtools depth -r chr18:0-80373285 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr18_depth_coverage_samtools.txt
samtools depth -r chr19:0-58617616 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr19_depth_coverage_samtools.txt
samtools depth -r chr20:0-64444167 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr20_depth_coverage_samtools.txt
samtools depth -r chr21:0-46709983 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr21_depth_coverage_samtools.txt
samtools depth -r chr22:0-50818468 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_Chr22_depth_coverage_samtools.txt
samtools depth -r chrX:0-156040895 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_ChrX_depth_coverage_samtools.txt
samtools depth -r chrY:0-57227415 --reference $hg38 $out/$input1/markdup/${nam1}_dedup.bam | awk '{sum+=$3} 
END{ print "Average = ",sum/NR}' > $out/$input1/alignments_stats/${input1}_ChrY_depth_coverage_samtools.txt

`grep "Average" $out/$input1/alignments_stats/${input1}_Chr*_depth_coverage_samtools.txt > $out/$input1/alignments_stats/All_Chr_depth_coverage_samtools.txt`


# Calculating Coverage for Tumor Sample for each chromosome

samtools depth -r chr1:0-248956422 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr1_depth_coverage_samtools.txt
samtools depth -r chr2:0-242193529 --reference $hg38  $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr2_depth_coverage_samtools.txt
samtools depth -r chr3:0-198295559 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr3_depth_coverage_samtools.txt
samtools depth -r chr4:0-190214555 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr4_depth_coverage_samtools.txt
samtools depth -r chr5:0-181538259 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr5_depth_coverage_samtools.txt
samtools depth -r chr6:0-170805979 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr6_depth_coverage_samtools.txt
samtools depth -r chr7:0-159345973 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr7_depth_coverage_samtools.txt
samtools depth -r chr8:0-145138636 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr8_depth_coverage_samtools.txt
samtools depth -r chr9:0-138394717 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr9_depth_coverage_samtools.txt
samtools depth -r chr10:0-133797422 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr10_depth_coverage_samtools.txt
samtools depth -r chr11:0-135086622 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr11_depth_coverage_samtools.txt
samtools depth -r chr12:0-133275309 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr12_depth_coverage_samtools.txt
samtools depth -r chr13:0-114364328 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr13_depth_coverage_samtools.txt
samtools depth -r chr14:0-107043718 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr14_depth_coverage_samtools.txt
samtools depth -r chr15:0-101991189 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3}END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr15_depth_coverage_samtools.txt
samtools depth -r chr16:0-90338345 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr16_depth_coverage_samtools.txt
samtools depth -r chr17:0-83257441 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr17_depth_coverage_samtools.txt
samtools depth -r chr18:0-80373285 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr18_depth_coverage_samtools.txt
samtools depth -r chr19:0-58617616 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr19_depth_coverage_samtools.txt
samtools depth -r chr20:0-64444167 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr20_depth_coverage_samtools.txt
samtools depth -r chr21:0-46709983 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr21_depth_coverage_samtools.txt
samtools depth -r chr22:0-50818468 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_Chr22_depth_coverage_samtools.txt
samtools depth -r chrX:0-156040895 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_ChrX_depth_coverage_samtools.txt
samtools depth -r chrY:0-57227415 --reference $hg38 $out/$input2/markdup/${nam2}_dedup.bam | awk '{sum+=$3} 
END{ print "Average = ",sum/NR}' > $out/$input2/alignments_stats/${input2}_ChrY_depth_coverage_samtools.txt

`grep "Average" $out/$input2/alignments_stats/${input2}_Chr*_depth_coverage_samtools.txt > $out/$input2/alignments_stats/All_Chr_depth_coverage_samtools.txt`


echo "Sambamba Alignment Statistics for $input1 and $input2 Completed" >> $out/WGS_analysis.log # Writing in log file
date >> $out/WGS_analysis.log									# Writing in log file
echo "##############################" >> $out/WGS_analysis.log					# Writing in log file

### Base Quality Score Recalibration For Normal

echo "gatk BaseRecalibrator --input $out/$input1/markdup/${nam1}_dedup.bam --reference $hg38 --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.known_indels.vcf --output $out/$input1/bqsr/${nam1}_recal_before_data.table" >> $out/parallel_BaseRecalibrator # Writing file for parallel run

#Apply BQSR For Normal
	
echo "gatk ApplyBQSR -R $hg38 -I $out/$input1/markdup/${nam1}_dedup.bam -bqsr $out/$input1/bqsr/${nam1}_recal_before_data.table -O $out/$input1/bqsr/${nam1}_dedup.bqsr.bam 2> $out/$input1/bqsr/${nam1}.applybqsr.stderr" >> $out/parallel_ApplyBQSR										# Writing file for parallel run

### Base Quality Score Recalibration For Tumor

echo "gatk BaseRecalibrator --input $out/$input2/markdup/${nam2}_dedup.bam --reference $hg38 --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.known_indels.vcf --output $out/$input2/bqsr/${nam2}_recal_before_data.table" >> $out/parallel_BaseRecalibrator

# Apply BQSR For Tumor

echo "gatk ApplyBQSR -R $hg38 -I $out/$input2/markdup/${nam2}_dedup.bam -bqsr $out/$input2/bqsr/${nam2}_recal_before_data.table -O $out/$input2/bqsr/${nam2}_dedup.bqsr.bam 2> $out/$input2/bqsr/${nam2}.applybqsr.stderr" >> $out/parallel_ApplyBQSR

echo "BQSR for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file

parallel -j 2 < $out/parallel_BaseRecalibrator	# Running script for baseRecalibrator for both normal and tumor in 
						# parallel

parallel -j 2 < $out/parallel_ApplyBQSR		# Running script fo ApplyBQSR for normal and tumor in parallel

echo "BQSR for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file
echo "##############################" >> $out/WGS_analysis.log		# Writing in log file


# Samtools indexing
echo "Samtools indexing for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log								# Writing in log file

	echo "samtools index -@ 60 $out/$input1/bqsr/${nam1}_dedup.bqsr.bam" >> $out/parallel_samtools_indexing
	
	echo "samtools index -@ 60 $out/$input2/bqsr/${nam2}_dedup.bqsr.bam" >> $out/parallel_samtools_indexing

parallel -j 2 < $out/parallel_samtools_indexing		# Running script in parallel

echo "Samtools indexing for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log								# Writing in log file
echo "##############################" >> $out/WGS_analysis.log				# Writing in log file


####### Module VariantCaller

# Creating folders for variant caller
mkdir -p $out/variantCaller/Varscan
mkdir -p $out/variantCaller/Mutect2
mkdir -p $out/variantCaller/lofreq
mkdir -p $out/variantCaller/strelka2_Germline
mkdir -p $out/variantCaller/HaplotypeCaller
mkdir -p $out/variantCaller/platypus

echo "First $input1 Second $input2 Third $input3"

# samtools mpileup
echo "samtools mpileup -B -f $hg38 $out/$input1/bqsr/${nam1}_dedup.bqsr.bam > $out/variantCaller/Varscan/${nam1}_dedup_bqsr_mpileup.bam" >> $out/parallel_mpileup

echo "samtools mpileup -B -f $hg38 $out/$input2/bqsr/${nam2}_dedup.bqsr.bam > $out/variantCaller/Varscan/${nam2}_dedup_bqsr_mpileup.bam" >> $out/parallel_mpileup

echo "Samtools mpileup $input1 and $input2 Started" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log

parallel -j 2 < $out/parallel_mpileup	# Running mpileup script in parallel

echo "Samtools mpileup $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file
echo "##############################" >> $out/WGS_analysis.log			# Writing in log file


echo "Varscan Started" >> $out/WGS_analysis.log       # Writing in log file
date >> $out/WGS_analysis.log

varscan somatic $out/variantCaller/Varscan/${nam1}_dedup_bqsr_mpileup.bam $out/variantCaller/Varscan/${nam2}_dedup_bqsr_mpileup.bam $out/variantCaller/Varscan/VS.vcf --output-vcf 1  # Running Varscan

echo "Varscan Completed" >> $out/WGS_analysis.log		# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file


### Segregating individual chromosome vcf and Creating Rscript for Varscan vcf reformating

for ((i=1; i <= 22; i++))

do

	cat $out/variantCaller/Varscan/VS.vcf.snp|grep "#" > $out/variantCaller/Varscan/Varscan_chr${i}.vcf.snp

	cat $out/variantCaller/Varscan/VS.vcf.snp|grep "^chr${i}	" >> $out/variantCaller/Varscan/Varscan_chr${i}.vcf.snp
	
	cat $out/variantCaller/Varscan/VS.vcf.indel|grep "#" > $out/variantCaller/Varscan/Varscan_chr${i}.vcf.indel

	cat $out/variantCaller/Varscan/VS.vcf.indel|grep "^chr${i}	" >> $out/variantCaller/Varscan/Varscan_chr${i}.vcf.indel

	echo "Rscript varscanReformat.R $out/variantCaller/Varscan/Varscan_chr${i}.vcf.snp" >> $out/parallel_Rscript
	
	echo "Rscript varscanReformat.R $out/variantCaller/Varscan/Varscan_chr${i}.vcf.indel" >> $out/parallel_Rscript

done

cat $out/variantCaller/Varscan/VS.vcf.snp|grep "#" > $out/variantCaller/Varscan/Varscan_chrX.vcf.snp

cat $out/variantCaller/Varscan/VS.vcf.snp|grep "^chrX	" >> $out/variantCaller/Varscan/Varscan_chrX.vcf.snp

cat $out/variantCaller/Varscan/VS.vcf.indel|grep "#" > $out/variantCaller/Varscan/Varscan_chrX.vcf.indel

cat $out/variantCaller/Varscan/VS.vcf.indel|grep "^chrX	" >> $out/variantCaller/Varscan/Varscan_chrX.vcf.indel

echo "Rscript varscanReformat.R $out/variantCaller/Varscan/Varscan_chrX.vcf.snp" >> $out/parallel_Rscript

echo "Rscript varscanReformat.R $out/variantCaller/Varscan/Varscan_chrX.vcf.indel" >> $out/parallel_Rscript

cat $out/variantCaller/Varscan/VS.vcf.snp|grep "#" > $out/variantCaller/Varscan/Varscan_chrY.vcf.snp

cat $out/variantCaller/Varscan/VS.vcf.snp|grep "^chrY	" >> $out/variantCaller/Varscan/Varscan_chrY.vcf.snp

cat $out/variantCaller/Varscan/VS.vcf.indel|grep "#" > $out/variantCaller/Varscan/Varscan_chrY.vcf.indel

cat $out/variantCaller/Varscan/VS.vcf.indel|grep "^chrY	" >> $out/variantCaller/Varscan/Varscan_chrY.vcf.indel

echo "Rscript varscanReformat.R $out/variantCaller/Varscan/Varscan_chrY.vcf.snp" >> $out/parallel_Rscript

echo "Rscript varscanReformat.R $out/variantCaller/Varscan/Varscan_chrY.vcf.indel" >> $out/parallel_Rscript


#Mutect2

echo "Mutect2 Started" >> $out/WGS_analysis.log		# Writing in log file
date >> $out/WGS_analysis.log				# Writing in log file

SECONDS=0

eval "$(conda shell.bash hook)"

conda activate wgs_gatk4


### Mutect2 running for each chromosome

pids=
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

	do
		gatk --java-options "-Xmx5g" Mutect2 -R $hg38 -I $out/$input2/bqsr/${nam2}_dedup.bqsr.bam -I $out/$input1/bqsr/${nam1}_dedup.bqsr.bam --germline-resource $db/GRCH38.P14/gnomAD/af-only-gnomad.hg38.vcf.gz --panel-of-normals $db/GRCH38.P14/panel-of-normals/1000g_pon.hg38.vcf.gz -O $out/variantCaller/Mutect2/$chr.vcf --native-pair-hmm-threads 10 -L $chr &
		pids+=" $!"
	done;

	echo "This will wait until all are done"

	date
	wait $pids
	wait
	duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed. "

echo "Mutect2 First step Done"

gatk --java-options "-Xmx16g" MergeVcfs \
	I=$out/variantCaller/Mutect2/chr1.vcf \
	I=$out/variantCaller/Mutect2/chr2.vcf \
	I=$out/variantCaller/Mutect2/chr3.vcf \
	I=$out/variantCaller/Mutect2/chr4.vcf \
	I=$out/variantCaller/Mutect2/chr5.vcf \
	I=$out/variantCaller/Mutect2/chr6.vcf \
	I=$out/variantCaller/Mutect2/chr7.vcf \
	I=$out/variantCaller/Mutect2/chr8.vcf \
	I=$out/variantCaller/Mutect2/chr9.vcf \
	I=$out/variantCaller/Mutect2/chr10.vcf \
	I=$out/variantCaller/Mutect2/chr11.vcf \
	I=$out/variantCaller/Mutect2/chr12.vcf \
	I=$out/variantCaller/Mutect2/chr13.vcf \
	I=$out/variantCaller/Mutect2/chr14.vcf \
	I=$out/variantCaller/Mutect2/chr15.vcf \
	I=$out/variantCaller/Mutect2/chr16.vcf \
	I=$out/variantCaller/Mutect2/chr17.vcf \
	I=$out/variantCaller/Mutect2/chr18.vcf \
	I=$out/variantCaller/Mutect2/chr19.vcf \
	I=$out/variantCaller/Mutect2/chr20.vcf \
	I=$out/variantCaller/Mutect2/chr21.vcf \
	I=$out/variantCaller/Mutect2/chr22.vcf \
	I=$out/variantCaller/Mutect2/chrX.vcf \
	I=$out/variantCaller/Mutect2/chrY.vcf \
	O=$out/variantCaller/Mutect2/Mutect2.vcf.gz

duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed. "


echo "Mutect2 Completed" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log

pids_vcf_reformat=
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
chr21 chr22 chrX chrY

do
	Rscript mutect2Reformat.R $out/variantCaller/Mutect2/$chr.vcf &
done

pids_vcf_reformat+=" $!"
wait $pids_vcf_reformat
wait

echo "Vcf reformat done"
date
duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed. "

echo "VCF reformat Done"

echo "Mutect2 Reformating Completed" >> $out/WGS_analysis.log			# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file
echo "##############################" >> $out/WGS_analysis.log		# Writing in log file


# HaplotypeCaller

echo "HaplotypeCaller Started" >> $out/WGS_analysis.log		# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file

SECONDS=0

eval "$(conda shell.bash hook)"

conda activate wgs_gatk4

### HaplotypeCaller running for each chromosome
pids=
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
      
do
	gatk --java-options "-Xmx5g" HaplotypeCaller -R $hg38 -I $out/$input2/bqsr/${nam2}_dedup.bqsr.bam --dbsnp $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -O $out/variantCaller/HaplotypeCaller/$chr.vcf --native-pair-hmm-threads 10 -L $chr &
      
	pids+=" $!"
done;

echo "This will wait until all  are done"
date
wait $pids
wait
duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed. "

echo "First step Done"

gatk --java-options "-Xmx16g" MergeVcfs \
      I=$out/variantCaller/HaplotypeCaller/chr1.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr2.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr3.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr4.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr5.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr6.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr7.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr8.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr9.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr10.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr11.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr12.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr13.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr14.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr15.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr16.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr17.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr18.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr19.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr20.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr21.vcf \
      I=$out/variantCaller/HaplotypeCaller/chr22.vcf \
      I=$out/variantCaller/HaplotypeCaller/chrX.vcf \
      I=$out/variantCaller/HaplotypeCaller/chrY.vcf \
      O=$out/variantCaller/HaplotypeCaller/${nam2}_HC.vcf

duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."


echo "HaplotypeCaller Completed" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log

pids_vcf_reformat=
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
chr21 chr22 chrX chrY

do
        Rscript haplotypeReformat.R $out/variantCaller/HaplotypeCaller/$chr.vcf &
done

pids_vcf_reformat+=" $!"
wait $pids_vcf_reformat
wait

echo "Vcf reformat done"
date
duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed. "

echo "VCF reformat Done"

echo "HaplotypeCaller Completed" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log
echo "##############################" >> $out/WGS_analysis.log

conda activate rna_preprocessing	# Activating conda environment


# Strelka
echo "Strelka Started" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log

/home/genomics/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam $out/$input2/bqsr/${nam2}_dedup.bqsr.bam --bam $out/$input2/bqsr/${nam2}_dedup.bqsr.bam --referenceFasta $hg38 --runDir $out/variantCaller/strelka2_Germline		# Creating files for Strelka to process

cd $out/variantCaller/strelka2_Germline		# Changing directory

python runWorkflow.py -m local -j 60		# Running python script for strelka


for ((i=1; i <= 22; i++))

do

        zcat $out/variantCaller/strelka2_Germline/results/variants/variants.vcf.gz|grep "#" > $out/variantCaller/strelka2_Germline/results/variants/Strelka_chr${i}.vcf

        zcat $out/variantCaller/strelka2_Germline/results/variants/variants.vcf.gz|grep "^chr${i}	" >> $out/variantCaller/strelka2_Germline/results/variants/Strelka_chr${i}.vcf

        echo "Rscript strelkaReformat.R $out/variantCaller/strelka2_Germline/results/variants/Varscan_chr${i}.vcf" >> $out/parallel_Rscript

done

zcat $out/variantCaller/strelka2_Germline/results/variants/variants.vcf.gz|grep "#" > $out/variantCaller/strelka2_Germline/results/variants/Strelka_chrX.vcf

zcat $out/variantCaller/strelka2_Germline/results/variants/variants.vcf.gz|grep "^chrX   " >> $out/variantCaller/strelka2_Germline/results/variants/Strelka_chrX.vcf

echo "Rscript strelkaReformat.R $out/variantCaller/strelka2_Germline/results/variants/Strelka_chrX.vcf" >> $out/parallel_Rscript

zcat $out/variantCaller/strelka2_Germline/results/variants/variants.vcf.gz|grep "#" > $out/variantCaller/strelka2_Germline/results/variants/Strelka_chrY.vcf

zcat $out/variantCaller/strelka2_Germline/results/variants/variants.vcf.gz|grep "^chrY	" > $out/variantCaller/strelka2_Germline/results/variants/Strelka_chrY.vcf

echo "Rscript strelkaReformat.R $out/variantCaller/strelka2_Germline/results/variants/Strelka_chrY.vcf" >> $out/parallel_Rscript

echo "Strelka Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log				# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file

cd $workdir	# Going back to working directory

conda activate rna_preprocessing


# Running platypus
echo "Platypus Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log				# Writing in log file

platypus callVariants --nCPU 60  --bamFiles $out/$input2/bqsr/${nam2}_dedup.bqsr.bam --refFile $hg38 --output $out/variantCaller/platypus/${nam2}_platypus.vcf	# Command for platypus


for ((i=1; i <= 22; i++))

do

        cat $out/variantCaller/platypus/${nam2}_platypus.vcf|grep "#" > $out/variantCaller/platypus/platypus_chr${i}.vcf

        cat $out/variantCaller/platypus/${nam2}_platypus.vcf|grep "^chr${i}	" >> $out/variantCaller/platypus/platypus_chr${i}.vcf

        echo "Rscript platypusReformat.R $out/variantCaller/platypus/platypus_chr${i}.vcf" >> $out/parallel_Rscript

done

cat $out/variantCaller/platypus/${nam2}_platypus.vcf|grep "#" > $out/variantCaller/platypus/platypus_chrX.vcf

cat $out/variantCaller/platypus/${nam2}_platypus.vcf|grep "^chrX   " >> $out/variantCaller/platypus/platypus_chrX.vcf

echo "Rscript platypusReformat.R $out/variantCaller/platypus/platypus_chrX.vcf" >> $out/parallel_Rscript

cat $out/variantCaller/platypus/${nam2}_platypus.vcf|grep "#" > $out/variantCaller/platypus/platypus_chrY.vcf

cat $out/variantCaller/platypus/${nam2}_platypus.vcf|grep "^chrY  " > $out/variantCaller/platypus/platypus_chrY.vcf

echo "Rscript platypusReformat.R  $out/variantCaller/platypus/platypus_chrY.vcf" >> $out/parallel_Rscript

echo "Platypus Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log				# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file

conda activate base

# Lofreq run
echo "Lofreq Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log			# Writing in log file

lofreq somatic -n $out/$input1/bqsr/${nam1}_dedup.bqsr.bam -t $out/$input2/bqsr/${nam2}_dedup.bqsr.bam -f $db/GRCH38.P14/hg38.fa --threads 60 -o $out/variantCaller/lofreq/out_ -d $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz # Command														# for lofreq
echo "Lofreq Completed" >> $out/WGS_analysis.log        # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log  # Writing in log file

### Segregating individual chromosome vcf and Creating Rscript for lofreq vcf reformating

zcat $out/variantCaller/lofreq/out_tumor_*snvs.vcf.gz|grep "^chr" > $out/variantCaller/lofreq/out_combined.snvs.vcf

zcat $out/variantCaller/lofreq/out_somatic_*snvs.vcf.gz|grep "^chr" >> $out/variantCaller/lofreq/out_combined.snvs.vcf

cat $out/variantCaller/lofreq/out_combined.snvs.vcf|sort -n|uniq > $out/variantCaller/lofreq/out_combined_uniq.snvs.vcf

zcat $out/variantCaller/lofreq/out_tumor_*indels.vcf.gz|grep "^chr" > $out/variantCaller/lofreq/out_combined.indels.vcf

zcat $out/variantCaller/lofreq/out_somatic_*indels.vcf.gz|grep "^chr" >> $out/variantCaller/lofreq/out_combined.indels.vcf

cat $out/variantCaller/lofreq/out_combined.indels.vcf|sort -n|uniq > $out/variantCaller/lofreq/out_combined_uniq.indels.vcf

for ((i=1; i <= 22; i++))
do
	zcat $out/variantCaller/lofreq/out_tumor_stringent.snvs.vcf.gz|grep "#" > $out/variantCaller/lofreq/lofreq_chr${i}.snvs.vcf
	
	cat $out/variantCaller/lofreq/out_combined_uniq.snvs.vcf|grep "^chr${i}	" >> $out/variantCaller/lofreq/lofreq_chr${i}.snvs.vcf

	zcat $out/variantCaller/lofreq/out_somatic_final.indels.vcf.gz|grep "#" > $out/variantCaller/lofreq/lofreq_chr${i}.indels.vcf
	
	cat $out/variantCaller/lofreq/out_combined_uniq.indels.vcf|grep "^chr${i}	" >> $out/variantCaller/lofreq/lofreq_chr${i}.indels.vcf

	echo "Rscript lofreqReformat.R $out/variantCaller/lofreq/lofreq_chr${i}.snvs.vcf" >> $out/parallel_Rscript
	
	echo "Rscript lofreqReformat.R $out/variantCaller/lofreq/lofreq_chr${i}.indels.vcf" >> $out/parallel_Rscript

done

zcat $out/variantCaller/lofreq/out_tumor_stringent.snvs.vcf.gz|grep "#" > $out/variantCaller/lofreq/lofreq_chrX.snvs.vcf

cat $out/variantCaller/lofreq/out_combined_uniq.snvs.vcf|grep "^chrX	" >> $out/variantCaller/lofreq/lofreq_chrX.snvs.vcf

zcat $out/variantCaller/lofreq/out_somatic_final.indels.vcf.gz|grep "#" > $out/variantCaller/lofreq/lofreq_chrX.indels.vcf

cat $out/variantCaller/lofreq/out_combined_uniq.indels.vcf|grep "^chrX	" >> $out/variantCaller/lofreq/lofreq_chrX.indels.vcf

echo "Rscript lofreqReformat.R $out/variantCaller/lofreq/lofreq_chrX.snvs.vcf" >> $out/parallel_Rscript

echo "Rscript lofreqReformat.R $out/variantCaller/lofreq/lofreq_chrX.indels.vcf" >> $out/parallel_Rscript

zcat $out/variantCaller/lofreq/out_tumor_stringent.snvs.vcf.gz|grep "#" > $out/variantCaller/lofreq/lofreq_chrY.snvs.vcf

cat $out/variantCaller/lofreq/out_combined_uniq.snvs.vcf|grep "^chrY	" >> $out/variantCaller/lofreq/lofreq_chrY.snvs.vcf

zcat $out/variantCaller/lofreq/out_somatic_final.indels.vcf.gz|grep "#" > $out/variantCaller/lofreq/lofreq_chrY.indels.vcf

cat $out/variantCaller/lofreq/out_combined_uniq.indels.vcf|grep "^chrY	" >> $out/variantCaller/lofreq/lofreq_chrY.indels.vcf

echo "Rscript lofreqReformat.R $out/variantCaller/lofreq/lofreq_chrY.snvs.vcf" >> $out/parallel_Rscript

echo "Rscript lofreqReformat.R $out/variantCaller/lofreq/lofreq_chrY.indels.vcf" >> $out/parallel_Rscript


echo "Lofreq Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log				# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file

conda activate wgs_gatk4

echo "Platypus, Strelka, Varscan and Lofreq reformating Started" >> $out/WGS_analysis.log        # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

parallel -j 55 < $out/parallel_Rscript

echo "Platypus, Strelka, Varscan and Lofreq reformating Completed" >> $out/WGS_analysis.log        # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log
####### Module Annovar

# Annovar analysis for Germline
echo "File prepare for Germline variant"

echo "Germline Annovar Started" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log

mkdir $out/Germline_Variant

grep "^chr" $out/variantCaller/HaplotypeCaller/*rf.vcf|cut -f 1,2 |sort -n |uniq > $out/Germline_Variant/HaplotypeCaller_location	# Fetching chromosome location from HaplotypeCaller vcf

cat $out/variantCaller/strelka2_Germline/results/variants/*rf.vcf|grep "^chr" |cut -f 1,2 |sort -n |uniq > $out/Germline_Variant/strelka2_germline_location	 # Fetching chromosome location from strelka vcf

grep "^chr" $out/variantCaller/platypus/*rf.vcf|cut -f 1,2 |sort -n |uniq > $out/Germline_Variant/platypus_location	#  Fetching chromosome location from platypus vcf

cat $out/Germline_Variant/*_location |sort -n |uniq -d -c |sort -r |grep -e "^      3" -e "^      2" > $out/Germline_Variant/concensus-list1		# Creating concensus list

echo "Concensus List & VCF making Starting"

sed -i 's/^      2 //g' $out/Germline_Variant/concensus-list1

sed -i 's/^      3 //g' $out/Germline_Variant/concensus-list1

grep "^chr" $out/variantCaller/HaplotypeCaller/*rf.vcf > $out/Germline_Variant/concensus-vcf
grep "^chr" $out/variantCaller/strelka2_Germline/results/variants/*rf.vcf > $out/Germline_Variant/concensus-vcf
grep "^chr" $out/variantCaller/platypus/*rf.vcf > $out/Germline_Variant/concensus-vcf

cp *.pl $out/		# Copying perl scripts

echo "1_concensus-vcf-Germline-find.pl running"

cd $out/

perl 1_concensus-vcf-Germline-find.pl	# Running concensus perl script

mkdir $out/Germline_Annovar	# Creating directory

cp $out/Germline_Variant/Final-Germline-ALL-CONCENSUS.vcf $out/Germline_Annovar/ 	# Copying Final concensus

echo "Starting of Germline Annovar"

cd $out/Germline_Annovar/
# Running Annovar
perl /home/genomics/software/annovar/table_annovar.pl Final-Germline-ALL-CONCENSUS.vcf $db/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 -out Variant-annovar.annotate -remove -protocol refGene,cytoBand,dbnsfp42a -operation g,r,f -nastring . -vcfinput

perl /home/genomics/software/annovar/annotate_variation.pl --geneanno Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype cosmic96_coding -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/Cosmic96_GRCH38/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype avsnp150 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype clinvar_20220320 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/table_annovar.pl -protocol dbnsfp42a -operation f -buildver hg38 Variant-annovar.annotate.avinput -out Variant-annovar.annotate.avinput.hg38_dbnsfp42a $db/annovar-database/humandb_July2022/GRCH38/humandb/ -nastring . -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype gme -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype exac03 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype tmcsnpdb -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype ljb26_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

sed -i 's/\t/\t#/' Variant-annovar.annotate.avinput.exonic_variant_function

grep -v "#synonymous" Variant-annovar.annotate.avinput.exonic_variant_function > Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

sed -i 's/\t#/\t/'  Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

echo "END of Germline Annovar"
cd $workdir

mkdir $out/RESULT_Germline
cd $out

echo "Tagging Germline Results"
perl 2_tagging_NewVCFFormat_Somatic.pl > $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info

sed -i 's/:/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/GT=//' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/;DP=/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/;RD=/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/;AD=/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/;AF=/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
cd $workdir

echo "Germline Annovar Completed" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log
echo "##############################" >> $out/WGS_analysis.log

################# Somatic variation
# Same for Somatic Annovar run
echo "Somatic Annovar Started" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log

mkdir $out/Somatic_Variant

grep "^chr" $out/variantCaller/Varscan/VS.vcf.* |cut -f 1,2 |sort -n |uniq > $out/Somatic_Variant/Varscan_location # Extracting genomic location from Varscan vcf

zcat $out/variantCaller/Mutect2/Mutect2.vcf.gz|grep "^chr" |cut -f 1,2 |sort -n |uniq > $out/Somatic_Variant/Mutect2_location	# Extracting genomic location from Mutect vcf

zcat $out/variantCaller/lofreq/out_somatic_final.*.vcf.gz|grep "^chr" | cut -f 1,2 |sort -n |uniq > $out/Somatic_Variant/lofreq_location	# Extracting genomic location from lofreq vcf

cat $out/Somatic_Variant/*_location |sort -n |uniq -d -c |sort -r |grep -e "^      3" -e "^      2" > $out/Somatic_Variant/concensus-list1

echo "Concensus List & VCF making Starting"

sed -i 's/^      2 //g' $out/Somatic_Variant/concensus-list1

sed -i 's/^      3 //g' $out/Somatic_Variant/concensus-list1

zcat $out/variantCaller/Mutect2/Mutect2.vcf.gz|grep "^chr" >> $out/Somatic_Variant/concensus-vcf

grep "^chr" $out/variantCaller/Varscan/VS.vcf* >> $out/Somatic_Variant/concensus-vcf

zcat $out/variantCaller/lofreq/out_somatic_final.*.vcf.gz|grep "^chr" >> $out/Somatic_Variant/concensus-vcf

echo "1_concensus-vcf-Germline-Somatic.pl running"

cd $out/

perl 1_concensus-vcf-Somatic-find.pl

mkdir $out/Somatic_Annovar

cp $out/Somatic_Variant/Final-Somatic-ALL-CONCENSUS.vcf $out/Somatic_Annovar/

echo "Starting of Somatic Annovar"

cd $out/Somatic_Annovar/
# Running Annovar for somatic
perl /home/genomics/software/annovar/table_annovar.pl Final-Somatic-ALL-CONCENSUS.vcf $db/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 -out Variant-annovar.annotate -remove -protocol refGene,cytoBand,dbnsfp42a -operation g,r,f -nastring . -vcfinput

perl /home/genomics/software/annovar/annotate_variation.pl --geneanno Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype cosmic96_coding -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/Cosmic96_GRCH38/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype avsnp150 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype clinvar_20220320 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/table_annovar.pl -protocol dbnsfp42a -operation f -buildver hg38 Variant-annovar.annotate.avinput -out Variant-annovar.annotate.avinput.hg38_dbnsfp42a $db/annovar-database/humandb_July2022/GRCH38/humandb/ -nastring . -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype gme -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype exac03 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype tmcsnpdb -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype ljb26_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo

sed -i 's/\t/\t#/' Variant-annovar.annotate.avinput.exonic_variant_function

grep -v "#synonymous" Variant-annovar.annotate.avinput.exonic_variant_function > Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

sed -i 's/#//' Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

echo "END of Somatic Annovar"

cd $workdir

mkdir $out/RESULT_Somatic

cd $out

echo "Tagging Results"

perl 2_tagging_NewVCFFormat_Somatic.pl > $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info 	# Writing files for 
											# somatic
sed -i 's/:/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/GT=//' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/;DP=/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/;RD=/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/;AD=/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/;AF=/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info

echo "Somatic Annovar Completed" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log
echo "##############################" >> $out/WGS_analysis.log


####### Module Copy number variation

# CNVpytor

conda activate base

echo "CNVpytor Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                   # Writing in log file

mkdir -p $out/cnv_output

cnvpytor -root $out/cnv_output/${sm2}_cnv.pytor -rd $out/$input2/bqsr/${nam2}_dedup.bqsr.bam -T $hg38

cnvpytor -root $out/cnv_output/${sm2}_cnv.pytor -his 10000 100000

cnvpytor -root $out/cnv_output/${sm2}_cnv.pytor -partition 10000 100000

cnvpytor -root $out/cnv_output/${sm2}_cnv.pytor -call 10000 > cnv_output/calls.10000.tsv

cnvpytor -root $out/cnv_output/${sm2}_cnv.pytor -call 100000 > cnv_output/calls.100000.tsv

cnvpytor -root $out/cnv_output/${sm2}_cnv.pytor -snp $out/variantCaller/HaplotypeCaller/${nam2}_HC.vcf -nofilter
#cnvpytor -root cnv_HC_output/file.pytor -pileup file.bam                   # OPTIONAL
#cnvpytor -root cnv_output/testSample.pytor -mask_snps                         # OPTIONAL
cnvpytor -root $out/cnv_output/${sm2}_cnv.pytor -baf 10000 100000

echo "CNVpytor Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                   # Writing in log file


####### Module Genefuse

# Genefuse on Tumor sample

echo "genefusion analysis started" >> $out/WGS_analysis.log
date >> $out/WGS_analysis.log

while IFS= read -r line
do

	echo $line
	header2=$(zcat $raw_tumor/$line | head -n 1);
        id2=$(echo $header2 | cut -f 3-4 -d":" | sed 's/@//');
        sm2=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");
        nam2=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-2 -d"_");


	if [ -z "$(ls -A  $out/$input2/trimmomatic_output/)" ]	# Checking folder if empty
	then

		genefuse -t 60 -r $hg38 -f $db/genefuse_database/cancer.hg38.csv -1 $raw_normal/${sm1}_R1.fastq -2 $raw_tumor/${sm2}_R2.fastq -u 5 -h $out/genefuse_output/res/report.html > $out/genefuse_output/res/result	# Run when trimmomatic_output is empty

	else
	
		genefuse -t 60 -r $hg38 -f $db/genefuse_database/cancer.hg38.csv -1 $out/$input1/trimmomatic_output/${sm1}_R1-trimmed_P.fastq -2 $out/$input2/trimmomatic_output/${sm2}_R2-trimmed_P.fastq -u 5 -h $out/genefuse_output/res/report.html > $out/genefuse_output/res/result	# Run when trimmomatic_output is not empty

	fi

	done < "$input4"

echo "genefusion analysis completed" >> $out/WGS_analysis.log		# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file


####### Module Microsatellite instability

# MSI Analysis

mkdir -p $out/MSI_output

echo "MSI analysis started" >> $out/WGS_analysis.log		# Writing in log file
date >> $out/WGS_analysis.log						                # Writing in log file

/home/genomics/software/msisensor/./msisensor -d /mnt/database/GRCH38.P14/hg38_homopolymer_microsatelittes.txt -n $out/$input1/bqsr/${nam1}_dedup.bqsr.bam -t $out/$input2/bqsr/${nam2}_dedup.bqsr.bam -b 60 -z 1 -o $out/MSI_output/${input1}_Vs_${input2}_msi.txt 

echo "MSI analysis completed" >> $out/WGS_analysis.log		# Writing in log file
date >> $out/WGS_analysis.log
echo "##############################" >> $out/WGS_analysis.log

####### Module Structural Variant

# Pindel Analysis
mkdir -p $out/Pindel_output

echo "Pindel 1st step analysis started" >> $out/WGS_analysis.log       # Writing in log file
date >> $out/WGS_analysis.log                                   # Writing in log file

samtools view -@ 60 $out/$input2/bqsr/${nam2}_dedup.bqsr.bam | /home/genomics/software/pindel/./sam2pindel - $out/Pindel_output/${nam2}_dedup.bqsr_Output4Pindel.txt 500 tumor 0 Illumina-PairEnd

echo "/home/genomics/software/pindel/pindel -f /mnt/database/GRCH38.P14/hg38.fa -p $out/Pindel_output/${nam2}_dedup.bqsr_Output4Pindel.txt -c ALL -T 60 -o $out/Pindel_output/${nam2}_dedup.bqsr_Final_output" >> $out/parallel_SV # Creating 
													#parallel_SV 
													#script for Pindel

echo "Pindel 1st step analysis completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file


# Delly Analysis

mkdir -p $out/Delly_output

echo "delly call -g $hg38 -o $out/Delly_output/delly.bcf $out/$input1/bqsr/${nam1}_dedup.bqsr.bam $out/$input2/bqsr/${nam2}_dedup.bqsr.bam" >> $out/parallel_SV		# Creating parallel_SV script for Delly

echo "Pindel and Delly analysis in parallel started" >> $out/WGS_analysis.log       # Writing in log file
date >> $out/WGS_analysis.log                                   # Writing in log file

parallel -j 2 < $out/parallel_SV	# Running parallel_SV; Pindel and Delly

echo "${nam2}_L001 tumor" >> $out/Delly_output/samples.tsv	# Writing sample file for Delly
echo "${nam1}_L001 control" >> $out/Delly_output/samples.tsv	# Writing sample file for Delly

delly filter -f somatic -o $out/Delly_output/delly.pre.bcf -s $out/Delly_output/samples.tsv $out/Delly_output/delly.pre.bcf/delly.bcf

conda deactivate

bcftools view $out/Delly_output/delly.pre.bcf > $out/Delly_output/delly_SV.vcf

grep "PASS" $out/Delly_output/delly_SV.vcf > $out/Delly_output/delly_SV_PASS.vcf

perl /home/genomics/software/annovar/table_annovar.pl delly_SV_PASS.vcf /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 -out $out/Delly_output/Variant-annovar.annotate -remove -protocol refGene,cytoBand,dbnsfp42a 
-operation g,r,f -nastring . -vcfinput

echo "Pindel and Delly analysis in parallel completed" >> $out/WGS_analysis.log       # Writing in log file
date >> $out/WGS_analysis.log                                   # Writing in log file
echo "##############################" >> $out/WGS_analysis.log  # Writing in log file

