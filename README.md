# SPPS Check

`spps_check.py` is a zero-dependency peptide feature extractor for quick human review and protein-design pipelines. It accepts raw sequences, FASTA files, PDB/mmCIF structures, directories, and stdin.

The concise views intentionally emphasize five first-pass features:

- peptide length
- longest hydrophobic run
- longest beta-branched run
- aspartimide-prone motif count
- cysteine count

These catch common coupling, on-resin aggregation, side-reaction, and disulfide-planning problems. They are not a guarantee of synthesis success. The JSON record also retains composition, physicochemical descriptors, all detected sites/runs, structure context, every triggered risk, and up to three ranked recommendations.

## Usage

```bash
# One sequence, human-readable
python spps_check.py KDIVPDDMQKAYELNKKRGENKLPTFEEVRSL

# One structure chain
python spps_check.py design.pdb --chain B

# Batch structure directory as one JSON object per line
python spps_check.py designs/ --kind structure --chain B --format jsonl > features.jsonl

# FASTA from stdin; fail a pipeline if any medium/high risk is found
cat designs.fasta | python spps_check.py - --format jsonl --fail-on medium

# Compact batch CSV/table
python spps_check.py designs.fasta --format csv
python spps_check.py designs/ --format table
```

Formats are `text`, `table`, `json`, `jsonl`, and `csv`. JSON formats contain the complete schema; table and CSV are deliberately compact. For mixed directories, use `--kind fasta` or `--kind structure`; `auto` intentionally processes every supported file. Exit codes are `0` for success, `1` for invalid input, and `2` when `--fail-on` is reached.

## Python API

```python
from spps_check import analyze_sequence, collect_records

record = analyze_sequence("ACDEFGHIK", name="candidate_1")
batch = collect_records(["designs/"])
```

`analyze_sequence` returns a JSON-serializable dictionary. `schema_version` changes if fields or meanings change incompatibly. Charge and pI assume free termini and unmodified standard amino acids; molecular weight is for the unmodified linear peptide.

## Default policy

Feature extraction is separate from risk policy. Thresholds are simple, inspectable heuristics in `_assess_risks`; projects should calibrate them against their resin, chemistry, peptide class, and acceptable failure rate. Recommendations are ranked by severity and capped at three, while the full `risks` array is never truncated.
