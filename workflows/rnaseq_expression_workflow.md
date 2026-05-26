1. Download mouse transcriptome from Ensembl
2. Build Salmon index:
   salmon index -t transcriptome/Mus_musculus.GRCm39.cdna.all.fa -i transcriptome/mouse_index

3. Download RNA-seq samples:
   prefetch SRR...

4. Convert SRA to FASTQ:
   fasterq-dump SRR.../SRR....sra --split-files

5. Compress FASTQ:
   gzip SRR..._1.fastq
   gzip SRR..._2.fastq

6. Quantify expression with Salmon:
   salmon quant -i transcriptome/mouse_index -l A -1 data/SRR..._1.fastq.gz -2 data/SRR..._2.fastq.gz -p 8 -o results/sample_quant

7. Extract selected gene TPMs:
   python scripts/extract_gene_tpm.py results/sample_quant/quant.sf results/sample_gene_tpm.tsv
