# pairwise_sim

A small, reproducible repository for simulating pairwise nucleotide alignment with **local BLAST+**.

The project is designed for simple experiments on one transcript or mRNA FASTA sequence. It supports:

- **sim1**: self-alignment of a sequence against itself
- **sim2**: alignment of a mutated sequence against the original sequence
- **sim3**: repeated alignment across a gradual error range, then downstream visualization in a notebook

The main alignment engine is **local BLAST+** (`blastn` + `makeblastdb`). Mutation generation is handled separately in Python so the intermediate mutated sequences can be inspected before alignment.

## Goals

- reproducible local workflow
- human-readable BLAST-style pairwise output
- saved intermediate FASTA files and mutation logs
- structured metrics for plotting and comparison
- easy to reuse with a different input FASTA

---

## Repository layout

```text
pairwise_sim/
├── README.md
├── environment.yml
├── .gitignore
├── data/
│   ├── raw/                # input FASTA files
│   └── generated/          # generated intermediate FASTA files
├── results/
│   ├── sim1/
│   ├── sim2/
│   └── sim3/
├── scripts/
│   ├── sim1.py
│   ├── sim2.py
│   ├── sim3.py
│   ├── mutate_sequence.py
│   ├── blast_runner.py
│   ├── parse_blast.py
│   └── utils.py
├── notebooks/
│   └── sim3_visualization.ipynb
├── docs/
│   └── example_outputs/
└── tests/
```

---

## Simulation definition

For this repository, mutation is defined as:

**error rate (%) = total edit events relative to the original sequence length**

For an original sequence of length `L`:

```text
n_edits = round(L * error_rate / 100)
```

Those edits are distributed across:

- substitutions
- insertions
- deletions

### Version 1 assumptions

- one substitution event changes **1 nucleotide**
- one insertion event inserts **1 nucleotide**
- one deletion event removes **1 nucleotide**
- the mutation fractions must sum to **1.0**
- mutation is applied relative to the original sequence
- random seed is user-controllable for reproducibility

This is a controlled simulation model for studying alignment behavior under noise. It is not intended to be a full biological evolution model.

---

## Input data

Put your input FASTA file in:

```text
data/raw/
```

Example expected usage:

```text
data/raw/NM_001178324.2.fa
```

You can replace that with any nucleotide FASTA file.

---

## Environment setup

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate pairwise-sim
```

The environment includes:

- Python 3.11
- BLAST+
- Biopython
- pandas
- matplotlib
- jupyterlab
- pytest

---

## Planned command-line interface

### sim1: self-alignment

```bash
python3 scripts/sim1.py \
  --input data/raw/NM_001178324.2.fa \
  --outdir results/sim1
```

### sim2: one error level

```bash
python3 scripts/sim2.py \
  --input data/raw/NM_001178324.2.fa \
  --error 20 \
  --sub 0.70 \
  --ins 0.15 \
  --del 0.15 \
  --seed 42 \
  --outdir results/sim2
```

### sim3: gradual error range

```bash
python3 scripts/sim3.py \
  --input data/raw/NM_001178324.2.fa \
  --min-error 1 \
  --max-error 20 \
  --step 1 \
  --sub 0.70 \
  --ins 0.15 \
  --del 0.15 \
  --seed 42 \
  --outdir results/sim3
```

---

## Expected outputs

### sim1

```text
results/sim1/<seqid>/
├── aln.txt
├── metrics.tsv
├── db/
└── logs/
```

### sim2

```text
results/sim2/<seqid>/e020_s42/
├── query_mut.fa
├── mut_log.tsv
├── aln.txt
├── metrics.tsv
├── db/
└── logs/
```

### sim3

```text
results/sim3/<seqid>/seed42/
├── runs/
│   ├── e001/
│   ├── e002/
│   └── ...
├── metrics_all.tsv
├── mut_summary.tsv
└── plots/
```

---

## Planned output files

### `query_mut.fa`
Mutated query FASTA generated before alignment.

### `mut_log.tsv`
Mutation event log.

Suggested columns:

```text
event_id	event_type	position_ref	ref_base	alt_base
```

Example:

```text
1	sub	104	A	G
2	del	287	T	-
3	ins	451	-	C
```

### `aln.txt`
Human-readable BLAST pairwise alignment report.

### `metrics.tsv`
Structured metrics table for downstream plotting.

Suggested columns:

```text
seq_id
error_rate
seed
len_ref
len_query
n_sub
n_ins
n_del
align_len
pct_ident
mismatch
gapopen
gaps
evalue
bitscore
score_raw
```

---

## Visualization ideas for sim3

The notebook will start with simple, interpretable plots such as:

- error rate vs bit score
- error rate vs percent identity
- error rate vs alignment length
- error rate vs number of gaps

Possible later additions:

- dot plot of reference vs mutated sequence
- heatmap for different mutation compositions
- repeated-run variability across multiple seeds

---

## Development plan

### Phase 1
- repository skeleton
- environment file
- README
- CLI placeholders

### Phase 2
- implement `sim1.py`
- create BLAST database locally
- run self-alignment and save pairwise output

### Phase 3
- implement `mutate_sequence.py`
- implement `sim2.py`
- save intermediate mutated FASTA and mutation log

### Phase 4
- implement `sim3.py`
- aggregate metrics across error rates
- build Jupyter notebook visualization

---

## Notes

- The default alignment program is expected to be `blastn` because the starting material is nucleotide sequence.
- BLAST+ supports local command-line workflows and custom database creation with `makeblastdb`.
- The repository is intentionally organized so another user can clone it, replace the FASTA file in `data/raw/`, and rerun the simulations.

## References

- NCBI BLAST+ command-line manual
- BLAST+ quick-start examples
- Bioconda `blast` package

-to be continued
