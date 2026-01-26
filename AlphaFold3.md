# AlphaFold3 on ALICE

---
## 1. input
---

AF3 accepts input in JSON format. For proteins, we have a script that converts fasta to json. It can take multiple fasta files at once if you refer to the directory containing the files. There are some default settings to the script that you can adjust once you understand what you are doing. 

Requirements for the script are
- Python3.10+
- Bio (python package, can be installed in conda using: conda install Bio -c bioconda)

```
python AF3_input_preparation.py --prey /path/to/directory/with/fasta/files --project name_of_project --make-sbatch --job-dirs

```
The full options of the script are a bit intimidating, but to begin you only need --prey, which refers to your fasta file, and --project, which is the directory into which your json file will go and your AlphaFold3 output later on. The --make-sbatch and --job-dirs options let the python script prepare the slurm job for you.
---
## 2. Submittig you slurm job
---

These are the historical default columns returned by `-outfmt 6`
when no custom field list is provided.

| #  | Field     | Description |
|----|-----------|-------------|
| 1  | qseqid    | Query sequence ID |
| 2  | sseqid    | Subject (target) sequence ID |
| 3  | pident    | Percentage of identical matches |
| 4  | length    | Alignment length |
| 5  | mismatch  | Number of mismatches |
| 6  | gapopen   | Number of gap openings |
| 7  | qstart    | Start of alignment in query |
| 8  | qend      | End of alignment in query |
| 9  | sstart    | Start of alignment in subject |
| 10 | send      | End of alignment in subject |
| 11 | evalue    | Expect value |
| 12 | bitscore  | Bit score |

This set is **not exhaustive** and is rarely optimal on its own.

---
## 3. Custom tabular output (recommended)
---
You can explicitly choose which columns BLAST writes:

    -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle"

---
## 4. Core blastp-specific fields
---

### Identifiers and annotation

- `qseqid`   : query sequence ID
- `sseqid`   : subject sequence ID
- `sacc`     : subject accession
- `saccver`  : subject accession.version
- `stitle`   : subject description line

### Alignment quality

- `pident`   : percent identity
- `length`   : alignment length
- `nident`   : number of identical residues
- `positive` : number of positive-scoring residues
- `gaps`     : total gaps
- `gapopen`  : gap openings

### Coordinates

- `qstart`, `qend`
- `sstart`, `send`

### Statistics

- `evalue`
- `bitscore`

---
## 5. Coverage and length fields (important)
---

Coverage fields are essential for filtering weak or spurious hits.

- `qlen`     : full query length
- `slen`     : full subject length
- `qcovs`    : query coverage per subject
- `qcovhsp`  : query coverage per HSP

Example:

    -outfmt "6 qseqid sacc pident length qlen qcovs evalue bitscore"

---
## 6. Taxonomy-aware fields (nr database)
---

These fields are especially useful when searching the NCBI nr database.

- `ssciname`   : subject scientific name
- `scomname`   : subject common name
- `staxids`    : NCBI taxonomy ID(s)
- `sskingdom`  : superkingdom

Example:

    -outfmt "6 qseqid sacc pident length qcovs evalue bitscore ssciname staxids"

---
## 7. Sequence-including fields (can infalte output)
---

These fields massively increase output size.

- `qseq`  : aligned query sequence
- `sseq`  : aligned subject sequence

Only use when manual alignment inspection is required.

---
## 8. Fields mainly for nucleotide-based BLAST
---

Relevant for blastn, tblastn, blastx, tblastx:

- `qframe`  : query frame
- `sframe`  : subject frame
- `sstrand` : subject strand

Irrelevant for blastp.

---
## 9. Recommended column sets
---

### General blastp annotation

    -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle"

### Strict homology filtering

    -outfmt "6 qseqid sacc pident length nident qlen qcovs evalue bitscore"

### Taxonomy-aware screening

    -outfmt "6 qseqid sacc pident length qcovs evalue bitscore ssciname staxids"

### Alignment inspection / debugging

    -outfmt "6 qseqid sacc qstart qend sstart send pident length qseq sseq"

---
## 10. Other BLAST programs (summary)
---

| Program  | Query type | Subject type | Notes |
|---------|------------|--------------|-------|
| blastp  | protein    | protein      | Most common for proteomics |
| blastn  | nucleotide | nucleotide   | Strand and frame fields relevant |
| blastx  | nucleotide | protein      | Query translated |
| tblastn | protein    | nucleotide   | Subject translated |
| tblastx | nucleotide | nucleotide   | Both translated |

Most fields are shared, but frame/strand fields only apply to nucleotide
or translated searches.

---
## 11. How to list all available fields locally
---

Run:

    blastp -help

Search for the section titled **Custom output formats**.
The exact field list may vary slightly between BLAST+ versions.

---
## 12. Final advice
---

- Always specify your `-outfmt` fields explicitly
- Always include at least one coverage metric
- Avoid sequence fields unless necessary
- Keep column order stable for downstream scripts

