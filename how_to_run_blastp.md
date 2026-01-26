# HPC ALICE BLASTP usage patterns

This guide shows common BLASTP use cases using built-in BLAST options.  

---

## Preparation before running BLAST

Activate the BLAST module
```
module load BLAST+/2.16.0-gompi-2024a
```
run this to ensure the database path is correctly used 
```
export BLASTDB=/zfsstore/databases/NCBI
```
To automate this for every login, add it to your bashrc

```
echo 'export BLASTDB=/zfsstore/databases/NCBI' >> ~/.bashrc
echo 'module load BLAST+/2.16.0-gompi-2024a' >> ~/.bashrc
source ~/.bashrc
```
---

## 1. Single protein vs entire nr database

Use this when you want to identify a protein or find homologs broadly.

    blastp \
      -query protein.faa \
      -db /zfsstore/databases/NCBI/nr/nr \
      -evalue 1e-5 \
      -max_target_seqs 20 \
      -num_threads 8 \
      -out protein_vs_nr.tsv \
      -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle"

---

## 2. Single protein vs a specific taxonomic group

Restrict results to a clade, order, genus, or species using NCBI taxids.  
NOTE: Taxonomy filtering limits biological scope, not database size.

    blastp \
      -query protein.faa \
      -db /zfsstore/databases/NCBI/nr/nr \
      -taxids 28211 \
      -evalue 1e-5 \
      -max_target_seqs 20 \
      -num_threads 8 \
      -out protein_vs_taxid.tsv \
      -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle"

---

## 3. Searching only a specific region of the protein

Use this when only part of the protein is biologically relevant.  
Residue numbering is 1-based and inclusive.

    blastp \
      -query protein.faa \
      -query_loc 120-240 \
      -db /zfsstore/databases/NCBI/nr/nr \
      -evalue 1e-5 \
      -max_target_seqs 20 \
      -num_threads 8 \
      -out protein_region_120_240.tsv \
      -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle"

---

## 4. Searching for divergent or remote homologs

Use this when close homologs are absent.  
This increases sensitivity at the cost of speed and specificity.

    blastp \
      -query protein.faa \
      -db /zfsstore/databases/NCBI/nr/nr \
      -word_size 2 \
      -matrix BLOSUM45 \
      -evalue 1e-2 \
      -qcov_hsp_perc 30 \
      -max_target_seqs 200 \
      -num_threads 8 \
      -out divergent_homologs.tsv \
      -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle"

---

## 5. Region-specific search for divergent homologs

Combines region restriction with relaxed similarity thresholds.

    blastp \
      -query protein.faa \
      -query_loc 120-240 \
      -db /zfsstore/databases/NCBI/nr/nr \
      -word_size 2 \
      -matrix BLOSUM45 \
      -evalue 1e-2 \
      -qcov_hsp_perc 30 \
      -max_target_seqs 200 \
      -num_threads 8 \
      -out region_divergent.tsv \
      -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle"

---

## 6. Performance and scaling guidelines

- Single protein: 1 job, 4–8 threads, <=15 minutes walltime
- Multiple proteins: use job arrays (one query per job)
- When setting cpus for your slurm job, make sure they are the same as number of threads used in blastp
- Avoid running BLAST on login nodes

---

## Summary of key options

Broad annotation        default blastp vs nr  
Taxonomy-restricted     -taxids  
Region-specific         -query_loc  
Divergent homologs      -word_size, -matrix, -evalue, -qcov_hsp_perc  
Throughput              job arrays + -num_threads
