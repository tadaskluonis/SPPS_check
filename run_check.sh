#!/usr/bin/env sh
python "$(dirname "$0")/spps_check.py" "$(dirname "$0")/pdb_and_fasta_files" --kind structure --chain B --format table
