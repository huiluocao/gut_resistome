#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mail-user=hcao@hku.hk
#SBATCH --mail-type=ALL
#SBATCH -p batch
#SBATCH --time=06-12:00:00

#### pls remind of the run time to be set

module use /sw/csi/modulefiles/applications;
module use /sw/csi/modulefiles/compilers;
module load ncbi-blast+/2.8.1/gnu-6.4.0;
module load bwa/0.7.17/gnu-6.4.0;
module load samtools/1.8 
module load bedtools/2.26.0/gnu-6.4.0;
module load perl;
module load python/2.7.14;

module load spades/3.11.1

##### input fastq file name and sample id
fastq_name=6014-17-001
sample_id=601417001

#### merge fastq files:
#zcat $fastq_name'_R1_001.fastq.gz' $fastq_name'_R1_002.fastq.gz' > $sample_id'_1.fastq.gz';
#zcat $fastq_name'_R2_001.fastq.gz' $fastq_name'_R2_002.fastq.gz' > $sample_id'_2.fastq.gz';

java -jar $TRIMMOMATIC_JAR PE -phred33 -summary $sample_id'_statsSummaryFile' -threads 16 $sample_id'_1.fastq.gz' $sample_id'_2.fastq.gz' \
$sample_id'_cleaned_1.fastq.gz' $sample_id'_unpaired_1.fastq.gz' \
$sample_id'_cleaned_2.fastq.gz' $sample_id'_unpaired_2.fastq.gz' \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

### bwa align against human host genome hg19 to remove host reads
#bwa mem /scratch/dragon/intel/x_caoh/microbiome_ref/pig_ref/GCA_000003025.6_Sscrofa11.1_genomic.fna $sample_id'_1.fastq.gz' $sample_id'_2.fastq.gz' > $sample_id.sam;
#bwa mem /scratch/dragon/intel/x_caoh/microbiome_ref/host_ref/hg19_rRNA_mito_Hsapiens_rna.fa $fastq_name'_1.fastq.gz' $fastq_name'_2.fastq.gz' > $sample_id.sam;
bwa mem /scratch/dragon/intel/x_caoh/microbiome_ref/broiler_ref/GCF_000002315.6_GRCg6a_genomic.fna \
raw_reads/$fastq_name'_1.fastq.gz' raw_reads/$fastq_name'_2.fastq.gz' > $sample_id.sam;
#### for pig, pls use '/scratch/dragon/intel/x_caoh/microbiome_ref/pig_ref/GCA_000003025.6_Sscrofa11.1_genomic.fna' as database
#### for chicken, pls use '/scratch/dragon/intel/x_caoh/microbiome_ref/broiler_ref/GCF_000002315.6_GRCg6a_genomic.fna' as reference
samtools view -S $sample_id.sam -b -o $sample_id.bam;
samtools view -b -f 4 $sample_id.bam > $sample_id'_unmapped.bam';
bamToFastq -i $sample_id'_unmapped.bam' -fq $sample_id'_unmapped_1.fastq' -fq2 $sample_id'_unmapped_2.fastq';
gzip $sample_id'_unmapped_1.fastq';
gzip $sample_id'_unmapped_2.fastq';
rm -rf $sample_id'_unmapped_1.fastq' $sample_id'_unmapped_2.fastq';
#rm -rf $sample_id'_1.fastq.gz' $sample_id'_2.fastq.gz';
rm -rf $sample_id.sam $sample_id.bam $sample_id'_unmapped.bam';
mv $sample_id'_unmapped_1.fastq.gz' cleaned_reads/;
mv $sample_id'_unmapped_2.fastq.gz' cleaned_reads/;

######assemble metagenome,sometime, if there is no enough mem in the node submitted, this step possibly can be failed. re-run is needed.
# calculate the total available RAM available for the allocated node
#mem=`sinfo -n $SLURMD_NODENAME -o %m --noheader`;
# convert memory from MB to GB
#mem=$((mem/1024));
spades.py --meta --debug --pe1-1 cleaned_reads/$sample_id'_unmapped_1.fastq.gz' --pe1-2 cleaned_reads/$sample_id'_unmapped_2.fastq.gz' \
-o ./$sample_id'_spades' -m 500 -k 21,33,55,77,85,99; ###--continue --only-assembler -k 21,33,55,77,99,127
perl /home/x_caoh/removesmalls.pl 100 ./$sample_id'_spades'/scaffolds.fasta > $sample_id'_contig100.fasta' ###remove short contigs <100bp
#### gene prediction and ARGs search
/scratch/dragon/intel/x_caoh/microbiome_ref/tools/MetaGeneMark_linux_64/mgm/gmhmmp $sample_id'_contig100.fasta' -r -m /scratch/dragon/intel/x_caoh/microbiome_ref/tools/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod -o $sample_id'.gff3' -f 3 -A $sample_id.faa -D $sample_id.ffn -V;

mv $sample_id'_contig100.fasta' contigs/;
mv $sample_id.ffn ffn/;
mv $sample_id.faa faa/;
mv $sample_id.gff3 gff3/;

bwa index ffn/$sample_id.ffn;
bwa mem ffn/$sample_id.ffn cleaned_reads/$sample_id'_unmapped_1.fastq.gz' cleaned_reads/$sample_id'_unmapped_2.fastq.gz' > $sample_id'_all_gene.sam';
#samtools view -S $sample_id'_all_gene.sam' -b -o $sample_id'_all_gene.bam';
#samtools view -b -F 4 $sample_id'_all_gene.bam' > $sample_id'_all_gene_mapped.bam';
/scratch/dragon/intel/x_caoh/microbiome_ref/tools/bbmap/pileup.sh in=$sample_id'_all_gene.sam' out=$sample_id'_coverage.txt' \
rpkm=$sample_id'_rpkm.txt' header=t headerpound=t covminscaf=10 minmapq=30;
mv $sample_id'_coverage.txt' ./coverage/;
mv $sample_id'_rpkm.txt' ./rpkm/;
rm -rf $sample_id'_all_gene.sam';

### new version of resfinder @ /scratch/dragon/intel/x_caoh/microbiome_ref/database/resfinder/2019May6_resfinder.ffn
blastn -query ffn/$sample_id.ffn -db /scratch/dragon/intel/x_caoh/microbiome_ref/database/resfinder/2019May6_resfinder.ffn \
-outfmt '6 qaccver saccver pident length positive ppos mismatch gapopen qstart qend sstart send qcovs evalue bitscore' \
-max_target_seqs 1 -out resfinder/$sample_id'_resfinder_args.txt' -evalue 1e-10;
#blastp -query faa/$sample_id.faa -db /scratch/dragon/intel/x_caoh/microbiome_ref/database/resfinder/resfinder_all_faa \
#-outfmt '6 qaccver saccver pident length positive ppos mismatch gapopen qstart qend sstart send qcovs evalue bitscore' \
#-max_target_seqs 1 -out resfinder/$sample_id'_resfinder_args.txt' -evalue 1e-10;
awk -F"\t" '{if (($6>=80)&&($13>=70)) {print$0}}' resfinder/$sample_id'_resfinder_args.txt' > resfinder/$sample_id'_arg_resfinder_filtered.txt';

#### metabat
module purge;
module load bwa/0.7.17/gnu-6.4.0;
module load samtools/1.8;
module load metabat/2.12.1.master/el7_gnu6.4.0_python2.7;
module load checkm/1.0.9/anaconda2-2.5.0;

bwa index ./contigs/$sample_id'_contig100.fasta';
bwa mem ./contigs/$sample_id'_contig100.fasta' cleaned_reads/$sample_id'_unmapped_1.fastq.gz' cleaned_reads/$sample_id'_unmapped_2.fastq.gz' > $sample_id'_contig100.sam';
samtools view -S $sample_id'_contig100.sam' -b -o $sample_id'_contig100.bam';
samtools sort $sample_id'_contig100.bam' -o $sample_id'_contig100_sorted.bam';
jgi_summarize_bam_contig_depths --outputDepth ./contigs/$sample_id'_contig100_depth.txt' $sample_id'_contig100_sorted.bam';
metabat2 -i ./contigs/$sample_id'_contig100.fasta' -a ./contigs/$sample_id'_contig100_depth.txt' -o metabin/$sample_id'_bins'/bin;
checkm lineage_wf -f metabin/$sample_id'_bins'/CheckM.txt -t 8 -x fa metabin/$sample_id'_bins'/ metabin/$sample_id'_bins'/SCG;
rm -rf $sample_id'_contig100.sam' $sample_id'_contig100.bam' $sample_id'_contig100_sorted.bam';

#### plasmid prediction using plasflow
module purge;
module load plasflow/1.1.0/anaconda2env
module load perl/5.22.4_gnu-640 ### to load bioperl

perl /scratch/dragon/intel/x_caoh/microbiome_ref/metaARG/China_chicken/scripts/filter_sequences_by_length.pl -input ./contigs/$sample_id'_contig100.fasta' \
-output ./contigs/$sample_id'_contig100_500more.fasta' -thresh 500;
PlasFlow.py --input ./contigs/$sample_id'_contig100_500more.fasta' --output plasmids/$sample_id'_plasflow_500_filtered';
cut -f1-6 plasmids/$sample_id'_plasflow_500_filtered' > plasmids/$sample_id'_plasflow_500_filtered_clean.txt';

#### plasmid prediction using recycler
module purge;
module load bwa/0.7.17/gnu-6.4.0;
module load samtools/1.8;
module load recycler/0.7/anaconda2env;

make_fasta_from_fastg.py -g $sample_id'_spades'/assembly_graph.fastg -o $sample_id'_spades'/assembly_graph.nodes.fasta;
bwa index $sample_id'_spades'/assembly_graph.nodes.fasta;
bwa mem $sample_id'_spades'/assembly_graph.nodes.fasta cleaned_reads/$sample_id'_unmapped_1.fastq.gz' cleaned_reads/$sample_id'_unmapped_2.fastq.gz' | samtools view -buS - > $sample_id'_spades'/reads_pe.bam;
samtools view -bF 0x0800 $sample_id'_spades'/reads_pe.bam > $sample_id'_spades'/reads_pe_primary.bam;
samtools sort $sample_id'_spades'/reads_pe_primary.bam -o $sample_id'_spades'/reads_pe_primary.sort.bam;
samtools index $sample_id'_spades'/reads_pe_primary.sort.bam;
recycle.py -g $sample_id'_spades'/assembly_graph.fastg -k 99 -b $sample_id'_spades'/reads_pe_primary.sort.bam -o plasmids/$sample_id'_recycler';

rm -rf $sample_id'_spades';

#### conjugation
module purge;
module load hmmer/3.2.1;
module load ncbi-blast+/2.8.1/gnu-6.4.0;
/home/x_caoh/anaconda2/bin/macsyfinder --db-type unordered -o conj/$sample_id'_conj' --sequence-db faa/$sample_id.faa \
-d /scratch/dragon/intel/x_caoh/microbiome_ref/metaARG/MGE/macsyfinder/Conjugation/DEF \
-p /scratch/dragon/intel/x_caoh/microbiome_ref/metaARG/MGE/macsyfinder/Conjugation/HMM \
--profile-suffix .HMM --hmmer /sw/csi/hmmer/3.2.1/el7_gnu6.4.0/bin/hmmsearch all


/home/x_caoh/anaconda2/bin/macsyfinder --db-type unordered -o CP023642_conj --sequence-db CP023642.faa -d /scratch/dragon/intel/x_caoh/KEEPME/microbiome/MGE_scripts/macsyfinder/Conjugation/DEF -p /scratch/dragon/intel/x_caoh/KEEPME/microbiome/MGE_scripts/macsyfinder/Conjugation/HMM --profile-suffix .HMM --hmmer /sw/csi/hmmer/3.2.1/el7_gnu6.4.0/bin/hmmsearch all


#### ARGs extract
/home/x_caoh/anaconda3/bin/python scripts/arg_get_v1.py $sample_id

### kraken2
/scratch/dragon/intel/x_caoh/microbiome_ref/tools/kraken2/kraken2 --db /scratch/dragon/intel/x_caoh/microbiome_ref/tools/kraken2/standard ./contigs/$sample_id'_contig100.fasta' --output kraken2/$sample_id'_contigs_kraken.txt' --use-names;

### metaphlan
module purge;
module load anaconda2/2.5.0 bowtie2/2.3.3.1;
zcat cleaned_reads/$sample_id'_unmapped_1.fastq.gz' cleaned_reads/$sample_id'_unmapped_2.fastq.gz' --to-stdout | /home/x_caoh/metaphlan/metaphlan.py --bowtie2db /home/x_caoh/metaphlan/bowtie2db/mpa --bt2_ps very-sensitive --input_type multifastq --bowtie2out metaphlan/$sample_id'_metaphlan.bt2out' > metaphlan/$sample_id'_metaphlan.outfmt6.txt';

### parser resfinder results with rpkm relative abundance
