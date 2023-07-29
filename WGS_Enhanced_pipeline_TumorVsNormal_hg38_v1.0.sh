#########################################################################################################################
#	This shell script takes input as folder name of Normal and Tumor.						#
#	e.g WGS_pipeline_TumorVsNormal_hg38_v0.1.sh Blood Tumor 							#
#	Last Modified on 12/03/2023											#
#	Last modified: Annovar correction										#
#															#
#															#
#															#
#	Tools Used in this pipeline											#
#	1.  Fastqc													#
#	2.  Trimmomatic													#
#	3.  BWA	mich													#
#	4.  Elprep													#
#	5.  HaplotypeCaller (Germline)											#
#       6.  Manta (Structural Variant)											#
#	7.  Strelka Germline (Germline)											#
#	8.  Strelka Somatic (Somatic)											#
#	9.  Annovar (Variant Annotation)										#
#########################################################################################################################

# Taking user inputs and setting resource paths
input1="$1"	# Taking input of normal folder name as argument
input2="$2"	# Taking input of tumor folder name as argument
workdir=$(pwd)	# storing path of Working directory
raw_normal="$workdir/$input1"	# Assigning path of normal sample
raw_tumor="$workdir/$input2"	# Assigning path of tumor sample 
out="$workdir/Output_hg38_${input1}_Vs_${input2}"	# Assigning path of output folder
db="/mnt/database"		# Assigning database folder path
hg38="$db/GRCH38.P14/hg38.fa"	# Assigning path of human genome reference file
hg38_mich="$db/GRCH38.P14/bwa_ert_database_hg38/hg38.fa"

# Removing pre-exiting files if any

rm $out/WGS_analysis.log
rm $out/parallel_fastqc
rm $out/parallel_trimmomatic
rm $out/parallel_after_trimmomatic
rm $out/parallel_BWA
rm $out/parallel_elprep
rm $out/parallel_BBIndex

# Checking files and folders Name format

if ls -d $input1|grep -E '^(KH-[0-9]{2})([A-Z]{2})([0-9]{1,3})-([A-Z]{2})_(S[0-9]{1,4})' && ls -d $input2|grep -E '^(KH-[0-9]{2})([A-Z]{2})([0-9]{1,3})-([A-Z]{2})_(S[0-9]{1,4})' && ls $input1|grep -E '^(KH-[0-9]{2})([A-Z]{2})([0-9]{1,3})-([A-Z]{2})_(S[0-9]{1,4})_(L[0-4]{3})_(R[1-2]{1})' && ls $input2|grep -E '^(KH-[0-9]{2})([A-Z]{2})([0-9]{1,3})-([A-Z]{2})_(S[0-9]{1,4})_(L[0-4]{3})_(R[1-2]{1})'
then
        echo Folder and File Name are okay
else
        echo ERROR: Please check Folder and File name
        exit 1
fi

# Merging fastq files and Creating list file
mkdir -p $out	# Creating Output folder
echo "Merging fastq files for $input1 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file

cd $input1

fn=$(ls | head -1| cut -f 1-2 -d "_")
echo $fn
zcat *R1* > ${fn}_All_R1_001.fastq
zcat *R2* > ${fn}_All_R2_001.fastq

ls *_All_R1_001* > ../list_${input1}

echo "Generated merged file size" >> $out/WGS_analysis.log	# Writing in log file
ls -trlh *All*  >> $out/WGS_analysis.log	# Writing in log file
echo "Merging fastq files for $input1 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file

cd $workdir

input3="list_${input1}"

echo "Merging fastq files for $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file
cd $input2

fn=$(ls | head -1| cut -f 1-2 -d "_")
echo $fn
zcat *R1* > ${fn}_All_R1_001.fastq
zcat *R2* > ${fn}_All_R2_001.fastq

ls *_All_R1_001* > ../list_${input2}


echo "Generated merged file size" >> $out/WGS_analysis.log	# Writing in log file
ls -trlh *All*  >> $out/WGS_analysis.log	# Writing in log file
echo "Merging fastq files for $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file

cd $workdir

input4="list_${input2}"


#Dynamic allocation of cpus

t="$(nproc --all)"      # Fetching total number of cpus
tnp=$((`expr $t - 3`)) # Maximum number of cpus will be utilized
val=$((`expr $tnp / 2`)) # Arithmatic operation on Max cpu count
np=$((`printf "%.*f\n" 0 $val`)) # Number of shared cpu



echo "First $input1 Second $input2 Third $input3 Forth $input4" # Print inputs as a test

# Creating folders for the intermediate results
mkdir -p $out/$input1/alignments_stats
mkdir -p $out/$input2/alignments_stats
mkdir -p $out/$input1/bqsr
mkdir -p $out/$input2/bqsr
mkdir -p $out/$input1/FastQC_output
mkdir -p $out/$input2/FastQC_output
mkdir -p $out/$input1/trimmomatic_output
mkdir -p $out/$input2/trimmomatic_output
mkdir -p $out/$input1/FastQC_output_after-trimmomatic
mkdir -p $out/$input2/FastQC_output_after-trimmomatic
mkdir -p $out/$input1/multiqc
mkdir -p $out/$input2/multiqc
mkdir -p $out/variantCaller/HaplotypeCaller
mkdir -p $out/variantCaller/Strelka2_Germline_Normal
mkdir -p $out/variantCaller/Strelka2_Germline_Tumor
mkdir -p $out/variantCaller/Strelka2_Somatic
mkdir -p $out/structuralVariant/Manta

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

		echo "fastqc $raw_normal/${sm}_R1* $raw_normal/${sm}_R2* -t $np -o $out/$input1/FastQC_output/" >> $out/parallel_fastqc

# Adaptor Trimming and cleaning

		echo "/home/genomics/anaconda3/envs/wgs_gatk4/bin/java -Xmx64g -Xmx64g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/trimmomatic-0.39-2/trimmomatic.jar PE -threads $np  -phred33 $raw_normal/${sm}_R1* $raw_normal/${sm}_R2* $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R1-trimmed_UP.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_UP.fastq ILLUMINACLIP:$db/Trimmomatic_adaptors/adaptors:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50" >> $out/parallel_trimmomatic

# QC-rechecking after trimming

		echo "fastqc -t $np $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq -o $out/$input1/FastQC_output_after-trimmomatic/" >> $out/parallel_after_trimmomatic

	done < "$input3"


# Fastqc and Trimmomatic run for Tumor

while IFS= read -r line
do
	echo $line
	sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");
	nam=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-2 -d"_");
	echo $sm
	echo $nam

# QC checking
		echo "fastqc $raw_tumor/${sm}_R1* $raw_tumor/${sm}_R2* -t $np -o $out/$input2/FastQC_output/" >> $out/parallel_fastqc
# Trimming
	
		echo "/home/genomics/anaconda3/envs/wgs_gatk4/bin/java -Xmx64g -Xmx64g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/trimmomatic-0.39-2/trimmomatic.jar PE -threads $np -phred33 $raw_tumor/${sm}_R1* $raw_tumor/${sm}_R2* $out/$input2/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input2/trimmomatic_output/${sm}_R1-trimmed_UP.fastq $out/$input2/trimmomatic_output/${sm}_R2-trimmed_P.fastq $out/$input2/trimmomatic_output/${sm}_R2-trimmed_UP.fastq ILLUMINACLIP:$db/Trimmomatic_adaptors/adaptors:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50" >> $out/parallel_trimmomatic

# QC-rechecking-after-trimming

		echo "fastqc -t $np $out/$input2/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input2/trimmomatic_output/${sm}_R2-trimmed_P.fastq -o $out/$input2/FastQC_output_after-trimmomatic/" >> $out/parallel_after_trimmomatic

done < "$input4"

echo "Fastqc for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file

parallel -j 2 < $out/parallel_fastqc		# Running script in parallel

echo "Fastqc for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file
echo "##############################" >> $out/WGS_analysis.log	# Writing in log file

echo "Trimmomatic for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log					# Writing in log file

parallel -j 2 < $out/parallel_trimmomatic		# Running script in parallel

echo "Trimmomatic for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file
echo "##############################" >> $out/WGS_analysis.log		# Writing in log file

echo "Fastqc after trimmomatic for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file

parallel -j 2 < $out/parallel_after_trimmomatic		# Running script in parallel

echo "Fastqc after trimmomatic for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file
echo "##############################" >> $out/WGS_analysis.log			# Writing in log file


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
		echo "/home/genomics/software/ert/./bwa-mem2 mem -K 100000000 -t $tnp -M -R \"@RG\\tID:sim\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:sim\\tPI:200\" -Z $hg38_mich $raw_normal/${sm1}_R1.fastq $raw_normal/${sm1}_R2.fastq 2> $out/$input1/${nam1}.align.stderr | samtools sort -@ $np -o $out/$input1/${nam1}.sorted.bam" >> $out/parallel_BWA	# Writing file for parallel run	when trimmomatic_output is empty
	else
		echo "/home/genomics/software/ert/./bwa-mem2 mem -K 100000000 -t $tnp -M -R \"@RG\\tID:sim\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:sim\\tPI:200\" -Z $hg38_mich $out/$input1/trimmomatic_output/${sm1}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm1}_R2-trimmed_P.fastq 2> $out/$input1/${nam1}.align.stderr | samtools sort -@ $np -o $out/$input1/${nam1}.sorted.bam" >> $out/parallel_BWA	# Writing file for 
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
		echo "/home/genomics/software/ert/./bwa-mem2 mem -K 100000000 -t $tnp -M -R \"@RG\\tID:sim\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:sim\\tPI:200\" -Z $hg38_mich $raw_tumor/${sm2}_R1.fastq $raw_tumor/${sm2}_R2.fastq 2> $out/$input2/${nam2}.align.stderr | samtools sort -@ $np -o $out/$input2/${nam2}.sorted.bam" >> $out/parallel_BWA
	else

		echo "/home/genomics/software/ert/./bwa-mem2 mem -K 100000000 -t $tnp -M -R \"@RG\\tID:sim\\tPL:ILLUMINA\\tLB:TruSeq\\tSM:sim\\tPI:200\" -Z $hg38_mich $out/$input2/trimmomatic_output/${sm2}_R1-trimmed_P.fastq $out/$input2/trimmomatic_output/${sm2}_R2-trimmed_P.fastq 2> $out/$input2/${nam2}.align.stderr | samtools sort -@ $np -o $out/$input2/${nam2}.sorted.bam" >> $out/parallel_BWA
	fi

	done < "$input4"

echo "BWA for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file

parallel -j 1 < $out/parallel_BWA	# Running script in parallel

echo "BWA for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log						# Writing in log file
echo "##############################" >> $out/WGS_analysis.log		# Writing in log file


###MarkDuplicates and BQSR For Normal

echo "elprep sfm $out/$input1/${nam1}.sorted.bam $out/$input1/bqsr/${nam1}_dedup_bqsr.bam --mark-duplicates --mark-optical-duplicates $out/$input1/bqsr/${nam1}_output.metrics --sorting-order coordinate  --bqsr $out/$input1/bqsr/${nam1}_output.recal  --reference /mnt/database/GRCH38.P14/elprep_database/hg38.elfasta --known-sites /mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.elsites,/mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.known_indels.elsites --nr-of-threads $tnp" > $out/parallel_elprep

###Cleaning caches

echo "sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'" >> $out/parallel_elprep

###MarkDuplicates and BQSR For Tumor

echo "elprep sfm $out/$input2/${nam2}.sorted.bam $out/$input2/bqsr/${nam2}_dedup_bqsr.bam --mark-duplicates --mark-optical-duplicates $out/$input2/bqsr/${nam2}_output.metrics --sorting-order coordinate  --bqsr $out/$input2/bqsr/${nam2}_output.recal --reference /mnt/database/GRCH38.P14/elprep_database/hg38.elfasta --known-sites /mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.elsites,/mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.known_indels.elsites --nr-of-threads $tnp --reference-confidence NONE --haplotypecaller $out/variantCaller/HaplotypeCaller/${nam2}_Haplotypecaller.vcf.gz" >> $out/parallel_elprep


echo "Elprep for $input1 and $input2 Started" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file

parallel -j 1 < $out/parallel_elprep	# Running script in parallel

sudo sh -c 'echo 3 >  /proc/sys/vm/drop_caches'

echo "Elprep for $input1 and $input2 Completed" >> $out/WGS_analysis.log	# Writing in log file
date >> $out/WGS_analysis.log							# Writing in log file
echo "##############################" >> $out/WGS_analysis.log			#


echo "BuildBamIndex for $input1 and $input2 Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                                                   # Writing in log file

echo "gatk --java-options \"-Xmx32G\" BuildBamIndex -I $out/$input1/bqsr/${nam1}_dedup_bqsr.bam -O $out/$input1/bqsr/${nam1}_dedup_bqsr.bai --MAX_RECORDS_IN_RAM 10000000" > $out/parallel_BBIndex
echo "gatk --java-options \"-Xmx32G\" BuildBamIndex -I $out/$input2/bqsr/${nam2}_dedup_bqsr.bam -O $out/$input2/bqsr/${nam2}_dedup_bqsr.bai --MAX_RECORDS_IN_RAM 10000000" >> $out/parallel_BBIndex

parallel -j 2 < $out/parallel_BBIndex

echo "BuildBamIndex for $input1 and $input2 Completed" >> $out/WGS_analysis.log        # Writing in log file
date >> $out/WGS_analysis.log                                                   # Writing in log file
echo "##############################" >> $out/WGS_analysis.log                  #


echo "Output copying to temporary folder Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

gsutil -m cp -r $out gs://kh-ngs-replica/Abhishek_bwa_mich/Spot_VM_backup_temporary/

echo "Output copying to temporary folder Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log


#SV Caller
eval "$(conda shell.bash hook)"
conda activate rna_preprocessing

echo "Manta analysis Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                                                   # Writing in log file

/home/genomics/software/manta-1.6.0.centos6_x86_64/bin/configManta.py --normalBam $out/$input1/bqsr/${nam1}_dedup_bqsr.bam --tumorBam $out/$input2/bqsr/${nam2}_dedup_bqsr.bam --referenceFasta $hg38 --runDir $out/structuralVariant/Manta

python $out/structuralVariant/Manta/runWorkflow.py -j $tnp 

echo "Manta analysis Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                                                   # Writing in log file
echo "##############################" >> $out/WGS_analysis.log

#Somatic Variant Caller

echo "Strelka Somatic analysis Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

/home/genomics/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam $out/$input1/bqsr/${nam1}_dedup_bqsr.bam --tumorBam $out/$input2/bqsr/${nam2}_dedup_bqsr.bam --referenceFasta $hg38 --indelCandidates $out/structuralVariant/Manta/results/variants/candidateSmallIndels.vcf.gz --runDir $out/variantCaller/Strelka2_Somatic


python $out/variantCaller/Strelka2_Somatic/runWorkflow.py -m local -j $tnp

echo "Strelka Somatic analysis Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log

#Germline Variant Caller

echo "Strelka Germline analysis for Normal Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

/home/genomics/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam $out/$input1/bqsr/${nam1}_dedup_bqsr.bam --referenceFasta $hg38 --runDir $out/variantCaller/Strelka2_Germline_Normal

python $out/variantCaller/Strelka2_Germline_Normal/runWorkflow.py -m local -j $tnp

echo "Strelka Germline analysis for Normal Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log

echo "Strelka Germline analysis for Tumor Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

/home/genomics/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --bam $out/$input2/bqsr/${nam2}_dedup_bqsr.bam --referenceFasta $hg38 --runDir $out/variantCaller/Strelka2_Germline_Tumor

python $out/variantCaller/Strelka2_Germline_Tumor/runWorkflow.py -m local -j $tnp

echo "Strelka Germline analysis for Tumor Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log

eval "$(conda shell.bash hook)"
conda activate base

#########CNV
./cnv_integrated.sh $nam $out


# Germline Annovar

mkdir -p $out/Germline_Annovar/

cp *.pl $out/

echo "Germline Annovar analysis Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

zcat $out/variantCaller/Strelka2_Germline_Tumor/results/variants/variants.vcf.gz|grep -w "PASS" > $out/Germline_Annovar/Pass_variants.vcf

python strelka2_germline_reformatting.py $out/Germline_Annovar/Pass_variants.vcf

cd $out/Germline_Annovar/

perl /home/genomics/software/annovar/table_annovar.pl Pass_variants_rf.vcf $db/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 --thread $tnp  --maxgenethread $tnp -out Variant-annovar.annotate -remove -protocol refGene,cytoBand,dbnsfp42a -operation g,r,f -nastring . -vcfinput

perl /home/genomics/software/annovar/annotate_variation.pl --geneanno Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 -otherinfo --thread $tnp --maxgenethread $tnp  --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype cosmic96_coding -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/Cosmic96_GRCH38/ -otherinfo --thread $tnp --maxgenethread $tnp  --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype avsnp150 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp  --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp  --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype clinvar_20220320 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/table_annovar.pl -protocol dbnsfp42a -operation f -buildver hg38 Variant-annovar.annotate.avinput -out Variant-annovar.annotate.avinput.hg38_dbnsfp42a $db/annovar-database/humandb_July2022/GRCH38/humandb/ -nastring . -otherinfo --thread $tnp

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype gme -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype exac03 -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype tmcsnpdb -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype ljb26_all -buildver hg38 Variant-annovar.annotate.avinput $db/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k

sed -i 's/\t/\t#/' Variant-annovar.annotate.avinput.exonic_variant_function

grep -v "#synonymous" Variant-annovar.annotate.avinput.exonic_variant_function > Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

sed -i 's/\t#/\t/'  Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

echo "END of Germline Annovar"
cd $workdir

mkdir $out/RESULT_Germline
cd $out

echo "Tagging Germline Results"
perl 2_tagging_NewVCFFormat_Germline.pl > $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info

sed -i 's/:/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/GT=//' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/;DP=/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/;RD=/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/;AD=/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info
sed -i 's/;AF=/\t/' $out/RESULT_Germline/${nam1}_${nam2}_Germline_Tagged-info

echo "Germline Annovar analysis Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log

cd $workdir

mkdir -p $out/Somatic_Annovar/

echo "Somatic Annovar analysis Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

zcat $out/variantCaller/Strelka2_Somatic/results/variants/somatic.indels.vcf.gz|grep -w "PASS" > $out/Somatic_Annovar/Pass_somatic_indels.vcf

python strelka2_somatic_indels_reformatting.py $out/Somatic_Annovar/Pass_somatic_indels.vcf

zcat $out/variantCaller/Strelka2_Somatic/results/variants/somatic.snvs.vcf.gz|grep -w "PASS" > $out/Somatic_Annovar/Pass_somatic_snvs.vcf

python strelka2_somatic_snvs_reformatting.py $out/Somatic_Annovar/Pass_somatic_snvs.vcf

cat $out/Somatic_Annovar/Pass_somatic_indels_rf.vcf > $out/Somatic_Annovar/Pass_somatic_rf.vcf

cat $out/Somatic_Annovar/Pass_somatic_snvs_rf.vcf >> $out/Somatic_Annovar/Pass_somatic_rf.vcf

cd $out/Somatic_Annovar/

perl /home/genomics/software/annovar/table_annovar.pl Pass_somatic_rf.vcf  /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 --thread $tnp --maxgenethread $tnp -out Variant-annovar.annotate -remove -protocol refGene,cytoBand,dbnsfp42a -operation g,r,f -nastring . -vcfinput

perl /home/genomics/software/annovar/annotate_variation.pl --geneanno Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype cosmic96_coding -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/Cosmic96_GRCH38/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype avsnp150 -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype clinvar_20220320 -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/table_annovar.pl -protocol dbnsfp42a -operation f -buildver hg38 Variant-annovar.annotate.avinput -out Variant-annovar.annotate.avinput.hg38_dbnsfp42a /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -nastring . -otherinfo --thread $tnp 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype gme -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype exac03 -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype tmcsnpdb -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype ljb26_all -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

sed -i 's/\t/\t#/' Variant-annovar.annotate.avinput.exonic_variant_function

grep -v "#synonymous" Variant-annovar.annotate.avinput.exonic_variant_function > Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

sed -i 's/#//' Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

echo "END of Somatic Annovar"
cd $workdir

mkdir $out/RESULT_Somatic

cd $out

echo "Tagging Somatic Results"
perl 2_tagging_NewVCFFormat_Somatic.pl > $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info

sed -i 's/:/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/GT=//' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/;DP=/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/;RD=/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/;AD=/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info
sed -i 's/;AF=/\t/' $out/RESULT_Somatic/${nam1}_${nam2}_Somatic_Tagged-info

echo "Somatic Annovar analysis Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log

#################################################################

cd $workdir

mkdir -p $out/Blood_Germline_Annovar/

echo "Germline from Blood Annovar analysis Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

zcat $out/variantCaller/Strelka2_Germline_Normal/results/variants/variants.vcf.gz|grep -w "PASS" > $out/Blood_Germline_Annovar/Pass_Blood_Germline.vcf

python strelka2_germline_reformatting.py $out/Blood_Germline_Annovar/Pass_Blood_Germline.vcf

cd $out/Blood_Germline_Annovar/

perl /home/genomics/software/annovar/table_annovar.pl Pass_Blood_Germline_rf.vcf /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 --thread $tnp --maxgenethread $tnp -out Variant-annovar.annotate -remove -protocol refGene,cytoBand,dbnsfp42a -operation g,r,f -nastring . -vcfinput

perl /home/genomics/software/annovar/annotate_variation.pl --geneanno Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -buildver hg38 -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype cosmic96_coding -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/Cosmic96_GRCH38/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype avsnp150 -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype clinvar_20220320 -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/table_annovar.pl -protocol dbnsfp42a -operation f -buildver hg38 Variant-annovar.annotate.avinput -out Variant-annovar.annotate.avinput.hg38_dbnsfp42a /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -nastring . -otherinfo --thread $tnp 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype gme -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype exac03 -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype tmcsnpdb -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

perl /home/genomics/software/annovar/annotate_variation.pl -filter -dbtype ljb26_all -buildver hg38 Variant-annovar.annotate.avinput /mnt/database/annovar-database/humandb_July2022/GRCH38/humandb/ -otherinfo --thread $tnp --maxgenethread $tnp --batchsize 80m --genomebinsize 500k 

sed -i 's/\t/\t#/' Variant-annovar.annotate.avinput.exonic_variant_function

grep -v "#synonymous" Variant-annovar.annotate.avinput.exonic_variant_function > Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

sed -i 's/#//' Variant-annovar.annotate.avinput.exonic_variant_function_non-synon

echo "END of Somatic Annovar"
cd $workdir

mkdir $out/RESULT_Blood_Germline

cd $out

echo "Tagging Somatic Results"
perl 2_tagging_NewVCFFormat_Germline.pl > $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info

sed -i 's/:/\t/' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info
sed -i 's/:/\t/' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info
sed -i 's/GT=//' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info
sed -i 's/;DP=/\t/' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info
sed -i 's/;RD=/\t/' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info
sed -i 's/;AD=/\t/' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info
sed -i 's/;AF=/\t/' $out/RESULT_Blood_Germline/${nam1}_Germline_Tagged-info

echo "Germline from Blood Annovar analysis Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log

# Copy to Storage
echo "Output folder copying Started" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file

gsutil -m cp -r $out gs://kh-ngs-replica/Abhishek_bwa_mich/Spot_VM_backup/

echo "Output folder copying Completed" >> $out/WGS_analysis.log  # Writing in log file
date >> $out/WGS_analysis.log                           # Writing in log file
echo "##############################" >> $out/WGS_analysis.log






