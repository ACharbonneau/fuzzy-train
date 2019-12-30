
module load bedtools
mkdir RawData
cd RawData

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

### PLINK files from https://www.nature.com/articles/s41588-018-0081-4

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz
tar -xvzf 1000G_Phase3_plinkfiles.tgz

## The hm3 SNP list doesn't have locations. So, to compare the output from using the
## hm3 SNPS vs all the SNPs in the subsets, I need to make it into a proper bed file
## This builds a SNP bed file from the recommended snplist: w_hm3.snplist and the location files
for i in `seq 1 22`
do cat 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim
done >> allchr.bim

cut -f 1 w_hm3.snplist > names_w_hm3.snplist
grep -Fwf names_w_hm3.snplist allchr.bim | awk '{ print "chr" $1 "\t" $4 "\t" $4+1 "\t" $2}' > w_hm3.bed


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

comm -12 RawData/UTR_3_UCSC.bed_LDSR.sort RawData/Coding_UCSC.bed_LDSR.sort > both3UTR_Coding.txt

# Use the master SNP bedfile from the LDSR paper to assign locations to SNPS and make a bed file

awk 'NR==FNR{a[$0]++;next}a[$4]'  both3UTR_Coding.txt RawData/w_hm3.bed > both3UTR_Coding.bed

# Find the intersection of the SNP bedfile and the splice bedfile

bedtools intersect -wa -wb -a both3UTR_Coding.bed -b RawData/VASTDB_PSI_Hsa108_hg19.bed > both3UTR_Coding.splice

# Sort the intersection file

sort -k 8 both3UTR_Coding.splice > both3UTR_Coding.splice.sorted

# Join the intersection file to the original splice file that has tissue information (makes an extremely wide excel sheet)

join -1 8 -2 2 both3UTR_Coding.splice.sorted RawData/VASTDB_PSI_Hsa108_hg19.tab > both3UTR_Coding.splice.tsv

# Add in headers

sed  -i '1i SpiceID\tSNPchr\tSNPstart\tSNPend\tSNPID\tSpliceChr\tSpliceStart\tSpliceEnd\tSpliceGene\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX\tAdipose_b\tAdipose_b-Q\tAdipose_c\tAdipose_c-Q\tAdipose_d\tAdipose_d-Q\tAdrenal_b\tAdrenal_b-Q\tAdrenal_c\tAdrenal_c-Q\tAmnion\tAmnion-Q\tAstrocytes\tAstrocytes-Q\tBladder_a\tBladder_a-Q\tBone_marrow_a\tBone_marrow_a-Q\tBone_marrow_b\tBone_marrow_b-Q\tBone_marrow_c\tBone_marrow_c-Q\tBrain_Endoth\tBrain_Endoth-Q\tBreast_Epith_a\tBreast_Epith_a-Q\tBreast_a\tBreast_a-Q\tCL_293T\tCL_293T-Q\tCL_Gm12878\tCL_Gm12878-Q\tCL_HeLa\tCL_HeLa-Q\tCL_K562\tCL_K562-Q\tCL_LP1\tCL_LP1-Q\tCL_MB231\tCL_MB231-Q\tCL_MCF7\tCL_MCF7-Q\tCL_PNT2\tCL_PNT2-Q\tCerebellum_a\tCerebellum_a-Q\tCerebellum_c\tCerebellum_c-Q\tChorion\tChorion-Q\tColon_b\tColon_b-Q\tColon_sigmoid\tColon_sigmoid-Q\tColon_transverse\tColon_transverse-Q\tCortex\tCortex-Q\tDecidua\tDecidua-Q\tESC_H1_a\tESC_H1_a-Q\tESC_H1_b\tESC_H1_b-Q\tESC_H1_c\tESC_H1_c-Q\tESC_H1_d\tESC_H1_d-Q\tESC_H9_a\tESC_H9_a-Q\tESC_H9_b\tESC_H9_b-Q\tEmbr_Cortex_13_17wpc\tEmbr_Cortex_13_17wpc-Q\tEmbr_Forebrain_9_12wpc\tEmbr_Forebrain_9_12wpc-Q\tEmbr_Forebrain_St13_14\tEmbr_Forebrain_St13_14-Q\tEmbr_Forebrain_St17_20\tEmbr_Forebrain_St17_20-Q\tEmbr_Forebrain_St22_23\tEmbr_Forebrain_St22_23-Q\tEndomStromCells\tEndomStromCells-Q\tEndothCells\tEndothCells-Q\tEpithelialCells\tEpithelialCells-Q\tFibroblasts\tFibroblasts-Q\tFrontal_Gyrus_old\tFrontal_Gyrus_old-Q\tFrontal_Gyrus_young\tFrontal_Gyrus_young-Q\tGLS_Cells\tGLS_Cells-Q\tHFDPC\tHFDPC-Q\tHMEpC_a\tHMEpC_a-Q\tHeart_a\tHeart_a-Q\tHeart_b\tHeart_b-Q\tHeart_c\tHeart_c-Q\tKidney_b\tKidney_b-Q\tKidney_c\tKidney_c-Q\tKidney_d\tKidney_d-Q\tLiver_a\tLiver_a-Q\tLiver_b\tLiver_b-Q\tLiver_c\tLiver_c-Q\tLung_b\tLung_b-Q\tLung_e\tLung_e-Q\tLung_f\tLung_f-Q\tLymph_node_b\tLymph_node_b-Q\tLymph_node_c\tLymph_node_c-Q\tMNC\tMNC-Q\tMSC\tMSC-Q\tMelanocytes\tMelanocytes-Q\tMicroglia\tMicroglia-Q\tMuscle_b\tMuscle_b-Q\tMuscle_d\tMuscle_d-Q\tMuscle_e\tMuscle_e-Q\tNPC_a\tNPC_a-Q\tNPC_b\tNPC_b-Q\tNeuroblastoma\tNeuroblastoma-Q\tNeurons\tNeurons-Q\tOligodendrocytes\tOligodendrocytes-Q\tOvary_a\tOvary_a-Q\tOvary_b\tOvary_b-Q\tPlacenta_Epith\tPlacenta_Epith-Q\tPlacenta_a\tPlacenta_a-Q\tPlacenta_b\tPlacenta_b-Q\tPlacenta_c\tPlacenta_c-Q\tProstate_b\tProstate_b-Q\tProstate_c\tProstate_c-Q\tProstate_d\tProstate_d-Q\tRetina_a\tRetina_a-Q\tRetina_macular\tRetina_macular-Q\tRetina_peripheral\tRetina_peripheral-Q\tSkin\tSkin-Q\tSmall_intestine\tSmall_intestine-Q\tSpleen_a\tSpleen_a-Q\tSpleen_b\tSpleen_b-Q\tStomach_a\tStomach_a-Q\tStomach_b\tStomach_b-Q\tSup_Temporal_Gyrus\tSup_Temporal_Gyrus-Q\tTestis_a\tTestis_a-Q\tTestis_b\tTestis_b-Q\tTestis_c\tTestis_c-Q\tThymus_a\tThymus_a-Q\tThymus_b\tThymus_b-Q\tThyroid_b\tThyroid_b-Q\tThyroid_c\tThyroid_c-Q\tThyroid_d\tThyroid_d-Q\tWBC_MNC_b\tWBC_MNC_b-Q\tWBC_MNC_c\tWBC_MNC_c-Q\tWhole_Brain_b\tWhole_Brain_b-Q\tiPS_a\tiPS_a-Q\tiPS_b\tiPS_b-Q' both3UTR_Coding.splice.tsv


