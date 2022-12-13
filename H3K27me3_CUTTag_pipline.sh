##############################  H3K27me3  ############################## 


mkdir {00.CleanData,01.RawData,02.QC,03.Bowtie2,04.Filter,05.ChIPseqSpikeInFree,06.SEACR,07.deepTools,08.MACS2,09.CHIPseeker,10.pyGenomeTracks,script}

### step 1

fastqc -t 6 -o ./02.QC ./00.CleanData/10+MET-1_1.clean.fq.gz ./00.CleanData/10+MET-1_2.clean.fq.gz
fastqc -t 6 -o ./02.QC ./00.CleanData/10+MET-2_1.clean.fq.gz ./00.CleanData/10+MET-2_2.clean.fq.gz
fastqc -t 6 -o ./02.QC ./00.CleanData/10MM-1_1.clean.fq.gz ./00.CleanData/10MM-1_2.clean.fq.gz
fastqc -t 6 -o ./02.QC ./00.CleanData/10MM-2_1.clean.fq.gz ./00.CleanData/10MM-2_2.clean.fq.gz
fastqc -t 6 -o ./02.QC ./00.CleanData/CON-1_1.clean.fq.gz ./00.CleanData/CON-1_2.clean.fq.gz
fastqc -t 6 -o ./02.QC ./00.CleanData/CON-2_1.clean.fq.gz ./00.CleanData/CON-2_2.clean.fq.gz

cd ./02.QC
multiqc .

bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa  Homo_sapiens.GRCh38.dna.primary_assembly


### step 2 alignment

for i in 10+MET-1 10+MET-2 10MM-1 10MM-2 CON-1 CON-2
do
bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 6 \
 -x /mnt/raid61/Personal_data/chenli/REF/GRCh38.93/bowtie2_index/Homo_sapiens.GRCh38.dna.primary_assembly \
 -1 ../00.CleanData/${i}_1.clean.fq.gz   \
 -2 ../00.CleanData/${i}_2.clean.fq.gz | samtools sort -O bam -@ 10 -o - > ${i}.bam 2> ./${i}_bowtie2.txt
done




### step 3 fragment Length

## Extract the 9th column from the alignment sam file which is the fragment length

for i in 10+MET-1 10+MET-2 10MM-1 10MM-2 CON-1 CON-2
do
samtools view -F 0x04 ./03.Bowtie2/${i}.bam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > ./03.Bowtie2/${i}_fragmentLen.txt

done




### step 4 filter


for i in 10+MET-1 10+MET-2 10MM-1 10MM-2 CON-1 CON-2
do
minQualityScore=2
samtools view -bh -q $minQualityScore ./03.Bowtie2/${i}.bam  > ./04.Filter/${i}.qualityScore$minQualityScore.bam
done



for i in 10+MET-1 10+MET-2 10MM-1 10MM-2 CON-1 CON-2
do

## Filter and keep the mapped read pairs
samtools view -bS -F 0x04 ./04.Filter/${i}.qualityScore2.bam > ./04.Filter/${i}.qualityScore2.mapped.bam


## sort and index
samtools sort -o ./04.Filter/${i}.qualityScore2.mapped.sorted.bam ./04.Filter/${i}.qualityScore2.mapped.bam                                                     
samtools index ./04.Filter/${i}.qualityScore2.mapped.sorted.bam 


## convert into bed file format
bedtools bamtobed -i ./04.Filter/${i}.qualityScore2.mapped.bam -bedpe > ./04.Filter/${i}.qualityScore2.mapped.bed


## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' ./04.Filter/${i}.qualityScore2.mapped.bed > ./04.Filter/${i}.qualityScore2.mapped.clean.bed 


## Only extract the fragment related columns
cut -f 1,2,6 ./04.Filter/${i}.qualityScore2.mapped.clean.bed | sort -k1,1 -k2,2n -k3,3n  > ./04.Filter/${i}.qualityScore2.mapped.fragments.bed

done






### step 5 normalized size factor
#### normalized using ChIPseqSpikeInFree (R package) 
#### /mnt/raid63/Cooperation/Liguideng_Tcell/X101SC21115238-Z01-J001/script/ChIPseqSpikeInFree_normalized.Rmd
#### output H3K27me3_SF.txt

for i in 10+MET-1 10+MET-2 10MM-1 10MM-2 CON-1 CON-2
do

SF=`grep ${i} ./05.ChIPseqSpikeInFree/H3K27me3_SF.txt|awk  '{print $7}'`
chromSize="/mnt/raid61/Personal_data/chenli/REF/GRCh38.93/hg38.chrom.sizes"

genomeCoverageBed -bg -scale $SF -ibam ./04.Filter/${i}.qualityScore2.mapped.sorted.bam  -g $chromSize | LC_COLLATE=C sort -k1,1 -k2,2n > ./05.ChIPseqSpikeInFree/${i}.normalized.bedGraph
bedGraphToBigWig ./05.ChIPseqSpikeInFree/${i}.normalized.bedGraph $chromSize ./05.ChIPseqSpikeInFree/${i}.normalized.bw

done






### step 6 SEACR

# for i in 10+MET-1 10+MET-2 10MM-1 10MM-2 CON-1 CON-2
# do

# seacr="/mnt/raid61/Personal_data/chenli/soft/SEACR-1.3/SEACR_1.3.sh"

# Calls enriched regions in target data by selecting the top 1% of regions by AUC

# bash $seacr ./05.SpikeInCalibration/${i}.normalized.bedgraph \
#     0.01 non stringent ./06.SEACR/${i}_seacr_top0.01.peaks

# done





### step 6 MACS2

for i in 10+MET-1 10+MET-2 10MM-1 10MM-2 CON-1 CON-2
do
macs2 callpeak -t ./04.Filter/${i}.qualityScore2.mapped.sorted.bam -g hs -f BAMPE --keep-dup all -n ${i}_macs2 --outdir ./08.MACS2 --broad --broad-cutoff 0.1 2> ./08.MACS2/${i}_macs2Peak_summary.txt

done




### step 7 deepTools


cat ./08.MACS2/10+MET-1_macs2_peaks.bed ./08.MACS2/10mm-1_macs2_peaks.bed ./08.MACS2/Ctrl-1_macs2_peaks.bed > ./08.MACS2/merge_macs2_peaks.bed

computeMatrix reference-point -S ./05.ChIPseqSpikeInFree/CON-1.normalized.bw ./05.ChIPseqSpikeInFree/10MM-1.normalized.bw  ./05.ChIPseqSpikeInFree/10+MET-1.normalized.bw \
                              -R ./08.MACS2/merge_macs2_peaks.bed \
                              --missingDataAsZero \
                              --beforeRegionStartLength 3000 \
                              --afterRegionStartLength 3000 \
                              --referencePoint center \
                              --skipZeros -o ./07.deepTools/MACS2/merge_MACS2.mat.gz -p 8


plotHeatmap -m ./07.deepTools/MACS2/merge_MACS2.mat.gz -out ./07.deepTools/MACS2/merge_MACS2_heatmap.png --sortUsing sum   --colorMap Reds


plotProfile -m ./07.deepTools/MACS2/merge_MACS2.mat.gz \
     -out ./07.deepTools/MACS2/merge_MACS2_line.pdf \
     --plotType lines --colors '#808080' '#FF6347' blue --perGroup \
     --plotFileFormat pdf \
     --plotWidth 10\
     --plotHeight 10\
     --samplesLabel Ctrl 10mM 10mM+Met




### step 8 CHIPseeker
#### using CHIPseeker R package 
#### /mnt/raid63/Cooperation/Liguideng_Tcell/X101SC21115238-Z01-J001/script/CHIPseeker_annalysis.Rmd

for i in 10+MET-1 10+MET-2 10MM-1 10MM-2 CON-1 CON-2
do

grep -v G ./08.MACS2/${i}_macs2_peaks.bed | grep -v K |grep -v MT > ../09.CHIPseeker/${i}_macs2_peaks_CHIPseeker.bed

sed -i "s/^/chr&/g" ../09.CHIPseeker/${i}_macs2_peaks_CHIPseeker.bed

done













