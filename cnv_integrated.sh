#read -e -p "please provide the location where your output directories are: " dir_in
#read -e -p "please provide the list of having your desired sample names: " sample

#if [ ! -f "$sample" ]; then
#        echo "bam list file not found"
#        exit 1
#fi
eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate base
workdir=$(pwd)

#cat $sample | while IFS= read -r nam
#do

#echo $nam
#cd $workdir
#mkdir $nam
nam=$1
dir=$2
#nam=$(sed -e 's/KH-22//' "$nam1" |sed -e 's/-.*//' |sed -e 's/_.*//')
echo $nam
echo $dir
mkdir -p $dir/CNV_Pytor

touch $dir/${nam}/bqsr/*.bai

cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -rd $dir/${nam}/bqsr/*.bam
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -his 10000 100000
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -partition 10000 100000
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -call 10000 > $dir/CNV_Pytor/$nam.10000.tsv
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -call 100000 > $dir/CNV_Pytor/$nam.100000.tsv
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -snp $dir/variantCaller/HaplotypeCaller/*$nam*.vcf.gz -nofilter
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -pileup $dir/${nam}/bqsr/*.bam                   # OPTIONAL
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -mask_snps                         # OPTIONAL
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -baf 10000 100000


cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -call combined 10000 > $dir/CNV_Pytor/${nam}.combined.10000.tsv
cnvpytor -j 40 -root $dir/CNV_Pytor/${nam}.pytor -call combined 100000 > $dir/CNV_Pytor/${nam}.combined.100000.tsv

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do

cnvpytor -root $dir/CNV_Pytor/${nam}.pytor -j 40 -view 100000 <<ENDL
set style classic
set rd_use_mask
set file_titles ${nam}_${chr}_in_1lakh
set panels rd likelihood snp
set markersize 0.2
$chr
save $dir/CNV_Pytor/${nam}_$chr_in_1lakh.png
ENDL
done

mkdir $dir/CNV_Pytor/${nam}_images
mv $dir/CNV_Pytor/*.png $dir/CNV_Pytor/${nam}_images


cnvpytor -root $dir/CNV_Pytor/${nam}.pytor -view 100000 <<ENDL
#set callers combined_mosaic        # IMPORTANT, default caller is mean shift
set Q0_range -1 0.5        # filter calls with more than half not uniquely mapped reads
set p_range 0 0.0001       # filter non-confident calls 
set pN_range 0 0.5              # filter calls with more than 50% Ns in reference genome 
set size_range 50000 inf   # filter calls smaller than 50kbp 
set dG_range 50000 inf    # filter calls close to gaps in reference genome (<100kbp)
print calls                # printing calls on screen (tsv format)
...
...
set print_filename $dir/CNV_Pytor/${nam}_hard_filtered_calls_100k.tsv   # output filename (xlsx, tsv or vcf)
set annotate               # turn on annotation (optional - takes a lot of time)
print calls                # generate output file with filtered calls 
quit

ENDL

cnvpytor -root $dir/CNV_Pytor/${nam}.pytor -view 10000 <<ENDL
#set callers combined_mosaic        # IMPORTANT, default caller is mean shift
set Q0_range -1 0.5        # filter calls with more than half not uniquely mapped reads
set p_range 0 0.0001       # filter non-confident calls
set pN_range 0 0.5              # filter calls with more than 50% Ns in reference genome
set size_range 50000 inf   # filter calls smaller than 50kbp
set dG_range 50000 inf    # filter calls close to gaps in reference genome (<100kbp)
print calls                # printing calls on screen (tsv format)
...
...
set print_filename $dir/CNV_Pytor/${nam}_hard_filtered_calls_10k.tsv   # output filename (xlsx, tsv or vcf)
set annotate               # turn on annotation (optional - takes a lot of time)
print calls                # generate output file with filtered calls
quit

ENDL

cnvpytor -root $dir/CNV_Pytor/${nam}.pytor -view 100000 <<ENDL
	set output_filename $dir/CNV_Pytor/${nam}_manhattan_plot_overall_100k.png
	set panels rd likelihood
	set rd_use_mask
	manhattan
ENDL

eval "$(conda shell.bash hook)"        # Setting bash for conda environment
#conda activate base
conda activate R_lang                # Activating conda environment
Rscript ${workdir}/R_script.R $dir/CNV_Pytor/${nam}_hard_filtered_calls_10k.tsv ${workdir}/cyto.txt
Rscript ${workdir}/R_script.R $dir/CNV_Pytor/${nam}_hard_filtered_calls_100k.tsv ${workdir}/cyto.txt

#rm *bqsr*
eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate base
#cd ..


#done < "$sample"
