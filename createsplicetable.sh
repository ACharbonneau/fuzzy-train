
module load bedtools
mkdir ../RawData
cd ../RawData

## Download files and uncompress

### original annotation files
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_v1.1_ldscores.tgz
tar -xvzf 1000G_Phase3_baseline_v1.1_ldscores.tgz

### splice locations
wget http://vastdb.crg.eu/downloads/Human/VASTDB_PSI_Hsa108_hg19.tab.gz
gunzip VASTDB_PSI_Hsa108_hg19.tab.gz

### Recommended SNP list from https://www.nature.com/articles/s41588-018-0081-4

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bzip2 -d w_hm3.snplist.bz2

# Reformat splice sites to a bed file

awk '{print $3 "\t" $2}' VASTDB_PSI_Hsa108_hg19.tab | sed 's/:/\t/g' | sed 's/-/\t/'  | awk  '$3!=""' |  awk  '$4!=""' |
tail -n +2  | sort | uniq > VASTDB_PSI_Hsa108_hg19.bed

# Extract annotation files

rm *_LDSR.txt
for col in `seq 6 52` 
    do for baseline in `seq 1 22`
        do NAME=`zcat ../RawData/baseline_v1.1/baseline.${baseline}.annot.gz | head -1 | cut -f ${col}`
        zcat ../RawData/baseline_v1.1/baseline.${baseline}.annot.gz | cut -f 3,${col} | awk '{if($2==1)print $1}' >> ${NAME}_LDSR.txt
    done 
done

# Make sorted versions of both annotation files and the splice file

sort Coding_UCSC.bed_LDSR.txt > Coding_UCSC.bed_LDSR.sort
sort UTR_3_UCSC.bed_LDSR.txt > UTR_3_UCSC.bed_LDSR.sort
sort -k 2 ../RawData/VASTDB_PSI_Hsa108_hg19.tab > ../RawData/VASTDB_PSI_Hsa108_hg19.tab.sorted


# Move out to main directory

cd ..

# Compare the two sorted annotation files, keep only lines that appear in both

comm -12 ../RawData/UTR_3_UCSC.bed_LDSR.sort ../RawData/Coding_UCSC.bed_LDSR.sort > both3UTR_Coding.txt

# Use the master SNP bedfile from the LDSR paper to assign locations to SNPS and make a bed file

awk 'NR==FNR{a[$0]++;next}a[$4]'  both3UTR_Coding.txt ../../10_2_Base/RawData/w_hm3.bed > both3UTR_Coding.bed

# Find the intersection of the SNP bedfile and the splice bedfile

bedtools intersect -wa -wb -a both3UTR_Coding.bed -b ../RawData/VASTDB_PSI_Hsa108_hg19.bed > both3UTR_Coding.splice

# Sort the intersection file

sort -k 8 both3UTR_Coding.splice > both3UTR_Coding.splice.sorted

# Join the intersection file to the original splice file that has tissue information (makes an extremely wide excel sheet)

join -1 8 -2 2 both3UTR_Coding.splice.sorted ../RawData/VASTDB_PSI_Hsa108_hg19.tab > both3UTR_Coding.splice.tsv



