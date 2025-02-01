**TASK:**
Blood samples from patients with SARS-CoV-2 have been received for analysis, paired samples (sample_1_1, sample_1_2 - sample_10_2). The goal is to determine the genome sequence by incorporating the detected variations into the reference genome Wuhan-Hu-1, evaluate the coverage of the S protein gene, determine the strain via the Pangolin server, and construct a phylogenetic tree.


**COMMANDS:**
Read quality qontrol:
```Shell
fasqc sample*_*.fastq -o /1st_fastq/
```
1.40% N content in short (up to 35 bp) reads
2.By the end, all reads contain G
3.Overall low quality, especially for short reads + very high number of short reads
4.Contamination is present

Trimming:
```Shell
for i in {1..10}; do fastp -i sample${i}_1.fastq -I sample${i}_2.fastq -o sample${i}1_trimmed.fastq -O sample${i}2_trimmed.fastq  --detect-adapter-for-pe  --trim_from1 20 --trim_from2 20 --trim_tail1 20 --trim_tail2 20 --cut-min-quality 20 --trim_poly_G   --n_base_limit 5  --length_required 40 --max_len1 150 --max_len2 150; done
fasqc sample*_*_trimmed.fastq  -o /trimmed/report
multiqc .
```

Indexing the reference:
```Shell
for i in {1..10}; do bowtie2 -x ~/project/wuhan_bowtie -1 ~/project/covid_raw/trimmed/sample${i}_1_trrimmed.fastq -2 ~/project/covid_raw/sample${i}_2_trimmed.fastq  -a  --very-sensitive  --local -N 1 -L 10 --rg "SM:sample${i}" --rg-id sample${i}  > ~/project/sorted/2/sample_${i}.bam; done
```

Mapping + Adding @RG:
```Shell
bowtie2 -a  --very-sensitive  --local -N 1 -L 10 --rg "SM:sample${i}" --rg-id sample${i}
```

Indexing and sorting:
```Shell
for i in {1..10};  do samtools sort ~/project/sorted/2/sample_${i}.bam -o ~/project/sorted/2/sample_${i}_sorted.bam; done
for i in {1..10}; do samtools index ~/project/sorted/2/sample_${i}_sorted.bam; done
for i in {1..10}; do samtools view sample_${i}_sorted.bam -o sample_${i}_sorted.sam; done
```

Examining the mapping results:
```Shell
parallel -j3 echo -ne {} "\"\\t\"" ";" samtools stats {} "|" grep ^SN "|" grep -P \"reads mapped\" -A 1  -B 1 ::: $(ls sample_*_sorted.bam)
```

Find close organisms:
```Shell
samtools stats sample_*_sorted.bam | grep ^SN
```

Obtaining coordinates of the S-gene from the annotation:
```Shell
grep -w 'gene "S"' wuhan.gtf | awk '$3=="gene"' > gene_S.gtf
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $3, ".", $7}' gene_S.gtf > gene_S.bed
```

Coverage assessment - creating .bed for the S-gene region from .gtf annotation:
```Shell
for i in {1..10}; do bedtools coverage -a ~/project/sorted/gene_S.bed -b ~/project/sorted/2/sample_${i}_sorted.bam >> coverage_gene_S.txt; done
```
The resulting file coverage_gene_S.txt contains the following columns:
1 – Chromosome
2 – Start
3 – End
4 – Region name (gene, stop, start codon, CDS)
5 – Score ('.' if not available)
6 – Strand
7 – Number of reads covering the region
8 – Region length
9 – Number of bases covered in the region
10 – Coverage fraction (from 0 to 1)

Where column 7 is the number of reads covering the region. Sum them using:
```Shell
awk '{sum += $7} END {print sum}' coverage_gene_S.txt
```
This gives the S-protein coverage: 71 230

Finding genomic variations:
```Shell
parallel -j 3 freebayes -f ~/project/wuhan.fasta -C 10 -F 0.5 {} ">" {.}.vcf ::: $(ls sample_*_sorted.bam)
parallel bgzip {} ::: $(ls *.vcf)
parallel tabix {} ::: $(ls *.vcf.gz)
vcf-merge -R *.vcf.gz > merged.vcf
```

Obtaining consensus sequences:
```Shell
for i in {1..10}; do bcftools consensus -f ~/project/wuhan.fasta ~/project/sorted/2/sample_${i}_sorted.vcf.gz -o ~/project/sorted/2/consensus/sample_${i}_consensus.fasta; done
```

Translating nucleotides into amino acid sequences + Extracting the S-gene:
```Shell
samtools faidx |  bcftools consensus  | sed "1d" | python -c "from Bio.Seq import Seq; import sys; print(f'>sample_{i}_S_protein\\n' + str(Seq(''.join([x.strip() for x in sys.stdin])).translate()))"
```

Extracting the S-gene from the consensus sequences:
```Shell
for i in {1..10}; do samtools faidx ~/project/sorted/2/consensus/sample_${i}_consensus.fasta MN908947.3:21563-25384 > ~/project/sorted/2/consensus/sample_${i}_S_gene.fasta; done
```

Combining the S-protein sequences with Known S-proteins from other strains:
```Shell
cat ~/project/sorted/2/consensus/S_protein_sequences/sample_*_S_protein.fasta > ~/project/sorted/2/consensus/S_protein_sequences/your_S_proteins.fasta
cp ~/project/covid_surfprot.pep.fa ~/project/sorted/2/consensus/S_protein_sequences/
cat ~/project/sorted/2/consensus/S_protein_sequences/your_S_proteins.fasta ~/project/sorted/2/consensus/S_protein_sequences/covid_surfprot.pep.fa > ~/project/sorted/2/consensus/S_protein_sequences/all_S_proteins.fasta
```

Translating to amino acid sequences:
```Shell
makeblastdb -in ~/project/sorted/2/consensus/S_protein_sequences/covid_surfprot.pep.fa -dbtype prot -out known_S_proteins_db
for i in {1..10}; do blastp -query ~/project/sorted/2/consensus/S_protein_sequences/sample_${i}_S_protein.fasta  -db known_S_proteins_db  -out ~/project/sorted/2/consensus/S_protein_sequences/sample_${i}_blastp.out  -outfmt 6  -evalue 1e-5  -num_threads 4;done
```

The resulting files consist of 12 columns, from which homologs can be selected:
1 – Query sequence ID: the name of the sequence
2 – Subject sequence ID: the name of the sequence from the database
3 – Percentage identity – should be over 95%
4 – Alignment length: the number of amino acids in the alignment – should be roughly that of the S protein (over 1000 AA)
5 – Number of mismatched amino acids
6 – Number of gap openings: insertions or deletions
7 – Start position (query)
8 – End position (query)
9 – Start position (subject)
10 – End position (subject)
11 – E-value
12 – Bit score – a measure of alignment quality

Count the number of homologs with high percentage identity:
```Shell
cat sample_${i}_blastp.out | awk '{if($3 > 90) print $0 }' | wc -l
```

Multiple sequence alignment:
```Shell
linsi ~/project/sorted/2/consensus/S_protein_sequences/all_S_proteins.fasta > ~/project/sorted/2/consensus/S_protein_sequences/all_S_proteins_aligned.fasta
```

Identification of amino acid substitutions in the multiple sequence alignment, present in several samples (>2):
```Shell
alv –only-variable all_S_proteins_aligned.fasta
```

Phylogenetic tree:
```Shell
FastTree all_S_proteins_aligned.fasta > S_proteins_aligned.tree
```


**RESULTS:**
Trimming performed:
*[link](https://drive.google.com/file/d/1e7lG9g4twnREYO8nXEShPGszEyXnsLCR/view?usp=sharing)*

Obtaining coordinates of the S-gene from the annotation:
*[link](https://drive.google.com/file/d/1Emeike-Y3-FfoEX2W16cWKUAF-8Gat2S/view?usp=sharing)*

Based on the analysis from the Pangolin website, the following strain data were obtained:
*[link](https://drive.google.com/file/d/1dpG9EZ1bDP4R569iMsQSoTMkoiND_vWX/view?usp=sharing)*
sample 1 — B.1.1.529  — Omicron
sample 2 — BA.1—Omicron
sample 3 — BA.1 —Omicron
sample 4 — B.1.617.2 —Delta
sample 5 — AY.45
sample 6 — AY.45
sample 7 — B.1
sample 8 — B.1.1.372 
sample 9 — B.1
sample 10 — B.

Phylogenetic tree:
*[link](https://drive.google.com/file/d/1a3k0itY1DWAhFIU5TbcecHD1rl-NSDXZ/view?usp=sharing)*

Presence of contamination: 
Using the command kraken2 and subsequent analysis of the report on the Pavian website:
```Shell
kraken2 --db /home/prep00/kraken_db/ –threads 4  --paired sample_*_1.fastq sample_*_2.fastq --report  kraken.rep --output kraken_output.txt
```

Detection of amino acid substitutions in the S-protein: 
Selecting only the S protein:
```Shell
samtools faidx wuhan.fasta MN908947.3:21563-25381 > S_sequence.fasta
```

Then, record:
```Shell
bcftools isec -n =1 -r MN908947.3:21563-25384 sample_1_sorted.vcf.gz sample_2_sorted.vcf.gz sample_3_sorted.vcf.gz sample_4_sorted.vcf.gz sample_5_sorted.vcf.gz sample_6_sorted.vcf.gz sample_7_sorted.vcf.gz sample_8_sorted.vcf.gz sample_9_sorted.vcf.gz sample_10_sorted.vcf.gz > variants_in_S.txt
```

Generate a file with nucleotide substitutions:
```Shell
awk '{print $1, $2-1, $2}' OFS="\t" variants_in_S.txt > coords.bed
```
MN908947.3      22450   C       T       0000001000
MN908947.3      22917   T       G       0001000000
MN908947.3      23012   G       A       0000001000
MN908947.3      23664   C       T       0000001000
MN908947.3      23731   C       T       0000000100
MN908947.3      24745   C       T       0001000000
MN908947.3      24795   C       T       0000001000
