#!/usr/bin/env python3
"""Lightweight SPPS feature extraction for peptide design pipelines."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Iterator, TextIO


SCHEMA_VERSION = "1.0"
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA_ORDER)
HYDROPHOBIC = set("AVILMFWYC")
BETA_BRANCHED = set("VIT")
ACIDIC = set("DE")
BASIC = set("KRH")

# Display priorities only. JSON output contains every calculated feature.
KEY_FEATURES = (
    "length",
    "max_hydrophobic_run",
    "max_beta_branched_run",
    "aspartimide_motif_count",
    "cysteine_count",
)
MAX_RECOMMENDATIONS = 3

ASPARTIMIDE_MOTIFS = ("DG", "DN", "DS", "DT", "DD")
DEAMIDATION_MOTIFS = ("NG", "NS", "NN", "NT")
ACID_LABILE_MOTIFS = ("DP",)

HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

# Average residue masses; water is added once for a linear unmodified peptide.
RESIDUE_MASS = {
    "A": 71.0788, "C": 103.1388, "D": 115.0886, "E": 129.1155,
    "F": 147.1766, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "K": 128.1741, "L": 113.1594, "M": 131.1926, "N": 114.1038,
    "P": 97.1167, "Q": 128.1307, "R": 156.1875, "S": 87.0782,
    "T": 101.1051, "V": 99.1326, "W": 186.2132, "Y": 163.1760,
}

# Approximate pKa model with free termini. Useful for comparative design features.
PKA_POSITIVE = {"n_term": 8.0, "K": 10.5, "R": 12.5, "H": 6.0}
PKA_NEGATIVE = {"c_term": 3.1, "D": 3.9, "E": 4.1, "C": 8.3, "Y": 10.1}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K",
}

SEVERITY_RANK = {"low": 1, "medium": 2, "high": 3}
RISK_PRIORITY = {
    "long_peptide": 0,
    "hydrophobic_run": 1,
    "beta_branched_run": 2,
    "aspartimide": 3,
    "cysteine_burden": 4,
}
SUPPORTED_SUFFIXES = {".fa", ".faa", ".fasta", ".pdb", ".cif", ".mmcif"}


def validate_sequence(sequence: str) -> str:
    """Normalize a sequence and reject, rather than silently delete, bad symbols."""
    sequence = "".join(sequence.split()).upper()
    if not sequence:
        raise ValueError("sequence is empty")
    invalid = sorted(set(sequence) - AA_SET)
    if invalid:
        raise ValueError(f"unsupported residue(s): {', '.join(invalid)}")
    return sequence


def _positions(sequence: str, residues: set[str]) -> list[int]:
    return [index + 1 for index, residue in enumerate(sequence) if residue in residues]


def _motif_positions(sequence: str, motifs: Iterable[str]) -> dict[str, list[int]]:
    return {
        motif: [
            index + 1
            for index in range(len(sequence) - len(motif) + 1)
            if sequence[index:index + len(motif)] == motif
        ]
        for motif in motifs
    }


def _runs(sequence: str, residues: set[str], minimum: int = 2) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    start = 0
    while start < len(sequence):
        if sequence[start] not in residues:
            start += 1
            continue
        end = start + 1
        while end < len(sequence) and sequence[end] in residues:
            end += 1
        if end - start >= minimum:
            result.append({
                "start": start + 1,
                "end": end,
                "length": end - start,
                "sequence": sequence[start:end],
            })
        start = end
    return result


def _charge(sequence: str, ph: float) -> float:
    counts = Counter(sequence)
    positive = 1.0 / (1.0 + 10.0 ** (ph - PKA_POSITIVE["n_term"]))
    negative = 1.0 / (1.0 + 10.0 ** (PKA_NEGATIVE["c_term"] - ph))
    for residue in ("K", "R", "H"):
        positive += counts[residue] / (1.0 + 10.0 ** (ph - PKA_POSITIVE[residue]))
    for residue in ("D", "E", "C", "Y"):
        negative += counts[residue] / (1.0 + 10.0 ** (PKA_NEGATIVE[residue] - ph))
    return positive - negative


def _isoelectric_point(sequence: str) -> float:
    low, high = 0.0, 14.0
    for _ in range(60):
        middle = (low + high) / 2.0
        if _charge(sequence, middle) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def _salt_bridge_candidates(sequence: str) -> list[dict[str, Any]]:
    result = []
    for index in range(len(sequence) - 4):
        left, right = sequence[index], sequence[index + 4]
        if (left in ACIDIC and right in BASIC) or (left in BASIC and right in ACIDIC):
            result.append({
                "start": index + 1,
                "end": index + 5,
                "pair": f"{left}{index + 1}-{right}{index + 5}",
            })
    return result


def _risk(
    risk_id: str,
    severity: str,
    feature: str,
    value: Any,
    threshold: str,
    message: str,
    recommendation: str,
) -> dict[str, Any]:
    return {
        "id": risk_id,
        "severity": severity,
        "feature": feature,
        "value": value,
        "threshold": threshold,
        "message": message,
        "recommendation": recommendation,
    }


def _assess_risks(
    features: dict[str, Any],
    structure: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    risks: list[dict[str, Any]] = []
    length = features["length"]
    hydrophobic_run = features["max_hydrophobic_run"]
    beta_run = features["max_beta_branched_run"]
    asp_count = features["aspartimide_motif_count"]
    cys_count = features["cysteine_count"]

    if length >= 60:
        risks.append(_risk("long_peptide", "high", "length", length, ">= 60",
            "Long linear peptides have a substantially larger cumulative coupling burden.",
            "Consider segment condensation or native chemical ligation."))
    elif length >= 40:
        risks.append(_risk("long_peptide", "medium", "length", length, ">= 40",
            "Length raises the cumulative probability of deletion products.",
            "Plan enhanced coupling/capping and verify intermediates during optimization."))

    if hydrophobic_run >= 5:
        risks.append(_risk("hydrophobic_run", "high", "max_hydrophobic_run", hydrophobic_run, ">= 5",
            "A long hydrophobic run can drive on-resin aggregation and poor solvation.",
            "Break the run if function permits, or plan aggregation-disrupting building blocks and solvent conditions."))
    elif hydrophobic_run >= 3:
        risks.append(_risk("hydrophobic_run", "medium", "max_hydrophobic_run", hydrophobic_run, ">= 3",
            "A hydrophobic run may reduce resin solvation and coupling efficiency.",
            "Inspect this region first when choosing resin, solvent, and coupling strategy."))

    if beta_run >= 4:
        risks.append(_risk("beta_branched_run", "high", "max_beta_branched_run", beta_run, ">= 4",
            "Consecutive beta-branched residues create a difficult coupling segment.",
            "Use stronger or repeated coupling around this segment and monitor deletions."))
    elif beta_run >= 3:
        risks.append(_risk("beta_branched_run", "medium", "max_beta_branched_run", beta_run, ">= 3",
            "Consecutive beta-branched residues can couple slowly.",
            "Plan stronger or repeated coupling around this segment."))

    if asp_count:
        risks.append(_risk("aspartimide", "medium", "aspartimide_motif_count", asp_count, ">= 1",
            "Asp-X motifs can form aspartimide during repeated base treatment.",
            "Review Asp protection and backbone-protection options at the reported sites."))

    if cys_count >= 4:
        risks.append(_risk("cysteine_burden", "high", "cysteine_count", cys_count, ">= 4",
            "Multiple cysteines increase protection and disulfide-isomer complexity.",
            "Define an orthogonal protection and oxidative-folding strategy before synthesis."))
    elif cys_count >= 2:
        risks.append(_risk("cysteine_burden", "medium", "cysteine_count", cys_count, ">= 2",
            "Cysteines require explicit thiol and disulfide handling.",
            "Plan cysteine protection and the intended disulfide-forming conditions."))
    elif cys_count == 1:
        risks.append(_risk("cysteine_handling", "low", "cysteine_count", cys_count, "= 1",
            "A free cysteine can oxidize or form mixed disulfides.",
            "Control thiol redox state during cleavage, purification, and storage."))

    if features["arginine_count"] >= 4:
        risks.append(_risk("arginine_rich", "medium", "arginine_count",
            features["arginine_count"], ">= 4",
            "Arginine-rich sequences can be demanding to deprotect and purify.",
            "Review Arg protecting-group removal and scavenger conditions."))

    if features["oxidation_sensitive_count"]:
        risks.append(_risk("oxidation_sensitive", "low", "oxidation_sensitive_count",
            features["oxidation_sensitive_count"], ">= 1",
            "Met and Trp can oxidize during synthesis, cleavage, purification, or storage.",
            "Limit oxidants and review cleavage scavengers and storage conditions."))

    if features["acid_labile_motif_count"]:
        risks.append(_risk("acid_labile_dp", "medium", "acid_labile_motif_count",
            features["acid_labile_motif_count"], ">= 1",
            "Asp-Pro bonds may cleave under strongly acidic conditions.",
            "Minimize unnecessary acid exposure and inspect cleavage products at the reported sites."))

    if features["gravy"] > 0 and abs(features["net_charge_ph7_4"]) < 1:
        risks.append(_risk("neutral_hydrophobic", "medium", "gravy/net_charge_ph7_4",
            {"gravy": features["gravy"], "net_charge_ph7_4": features["net_charge_ph7_4"]},
            "GRAVY > 0 and |charge| < 1",
            "A neutral, hydrophobic peptide may aggregate and be difficult to purify.",
            "Review solubilization and purification conditions, or add charge if function permits."))

    if features["deamidation_motif_count"]:
        risks.append(_risk("deamidation", "low", "deamidation_motif_count",
            features["deamidation_motif_count"], ">= 1",
            "Asn-X motifs can deamidate during handling or storage.",
            "Control pH, temperature, and storage time for the final peptide."))

    if structure and structure.get("beta_sheet_records", 0):
        risks.append(_risk("beta_sheet_structure", "low", "beta_sheet_records",
            structure["beta_sheet_records"], ">= 1",
            "The supplied structure has beta-sheet annotations, suggesting possible post-cleavage aggregation.",
            "Assess final-peptide solubility and aggregation experimentally."))

    return sorted(risks, key=lambda item: (
        -SEVERITY_RANK[item["severity"]],
        RISK_PRIORITY.get(item["id"], 99),
        item["id"],
    ))


def analyze_sequence(
    sequence: str,
    name: str = "sequence",
    structure_path: str | None = None,
    *,
    source: str | None = None,
    chain: str | None = None,
    structure: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return one stable, JSON-serializable SPPS feature record."""
    sequence = validate_sequence(sequence)
    counts = Counter(sequence)
    length = len(sequence)
    hydrophobic_runs = _runs(sequence, HYDROPHOBIC)
    beta_runs = _runs(sequence, BETA_BRANCHED)
    asp_sites = _motif_positions(sequence, ASPARTIMIDE_MOTIFS)
    deamidation_sites = _motif_positions(sequence, DEAMIDATION_MOTIFS)
    acid_labile_sites = _motif_positions(sequence, ACID_LABILE_MOTIFS)
    pseudo_sites = [
        {"position": index + 1, "dipeptide": sequence[index:index + 2]}
        for index in range(length - 1)
        if sequence[index + 1] in {"S", "T"}
    ]

    hydrophobic_count = sum(counts[aa] for aa in HYDROPHOBIC)
    beta_count = sum(counts[aa] for aa in BETA_BRANCHED)
    charged_count = sum(counts[aa] for aa in ACIDIC | BASIC)
    sites = {
        "aspartimide_motifs": asp_sites,
        "deamidation_motifs": deamidation_sites,
        "acid_labile_motifs": acid_labile_sites,
        "cysteine": _positions(sequence, {"C"}),
        "methionine": _positions(sequence, {"M"}),
        "tryptophan": _positions(sequence, {"W"}),
        "arginine": _positions(sequence, {"R"}),
        "proline": _positions(sequence, {"P"}),
        "glycine": _positions(sequence, {"G"}),
        "pseudoproline_candidates": pseudo_sites,
        "i_i4_salt_bridge_candidates": _salt_bridge_candidates(sequence),
    }
    features = {
        "length": length,
        "molecular_weight_da": round(sum(RESIDUE_MASS[aa] for aa in sequence) + 18.01528, 3),
        "gravy": round(sum(HYDROPATHY[aa] for aa in sequence) / length, 4),
        "isoelectric_point": round(_isoelectric_point(sequence), 3),
        "net_charge_ph1": round(_charge(sequence, 1.0), 3),
        "net_charge_ph7_4": round(_charge(sequence, 7.4), 3),
        "net_charge_ph11": round(_charge(sequence, 11.0), 3),
        "hydrophobic_count": hydrophobic_count,
        "hydrophobic_fraction": round(hydrophobic_count / length, 4),
        "max_hydrophobic_run": max((run["length"] for run in hydrophobic_runs), default=0),
        "beta_branched_count": beta_count,
        "beta_branched_fraction": round(beta_count / length, 4),
        "max_beta_branched_run": max((run["length"] for run in beta_runs), default=0),
        "charged_count": charged_count,
        "charged_fraction": round(charged_count / length, 4),
        "acidic_count": sum(counts[aa] for aa in ACIDIC),
        "basic_count": sum(counts[aa] for aa in BASIC),
        "arginine_count": counts["R"],
        "cysteine_count": counts["C"],
        "methionine_count": counts["M"],
        "tryptophan_count": counts["W"],
        "oxidation_sensitive_count": counts["M"] + counts["W"],
        "proline_count": counts["P"],
        "glycine_count": counts["G"],
        "aspartimide_motif_count": sum(map(len, asp_sites.values())),
        "deamidation_motif_count": sum(map(len, deamidation_sites.values())),
        "acid_labile_motif_count": sum(map(len, acid_labile_sites.values())),
        "pseudoproline_candidate_count": len(pseudo_sites),
        "i_i4_salt_bridge_candidate_count": len(sites["i_i4_salt_bridge_candidates"]),
        "residue_counts": {aa: counts[aa] for aa in AA_ORDER},
        "residue_fractions": {aa: round(counts[aa] / length, 4) for aa in AA_ORDER},
    }

    # Compatibility for callers of the old positional API.
    if structure is None and structure_path:
        path = Path(structure_path)
        structure = {"format": path.suffix.lstrip("."), "beta_sheet_records": _pdb_sheet_count(path)}

    risks = _assess_risks(features, structure)
    recommendations: list[str] = []
    for item in risks:
        if item["recommendation"] not in recommendations:
            recommendations.append(item["recommendation"])
        if len(recommendations) == MAX_RECOMMENDATIONS:
            break
    if not recommendations:
        recommendations = ["No major sequence-level SPPS liabilities detected by the default policy."]

    return {
        "schema_version": SCHEMA_VERSION,
        "id": name,
        "source": source,
        "chain": chain,
        "sequence": sequence,
        "features": features,
        "key_features": {key: features[key] for key in KEY_FEATURES},
        "runs": {"hydrophobic": hydrophobic_runs, "beta_branched": beta_runs},
        "sites": sites,
        "structure": structure,
        "risk_level": risks[0]["severity"] if risks else "none",
        "risks": risks,
        "recommendations": recommendations,
    }


def _parse_fasta_text(text: str, source: str) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name = line[1:].strip() or f"record_{line_number}"
            chunks = []
        else:
            if name is None:
                raise ValueError(f"{source}:{line_number}: sequence found before FASTA header")
            chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def _pdb_sheet_count(path: Path) -> int:
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            return sum(1 for line in handle if line.startswith("SHEET "))
    except OSError:
        return 0


def _parse_pdb(path: Path, selected_chain: str | None) -> Iterator[tuple[str, str, dict[str, Any]]]:
    residues: dict[str, list[str]] = {}
    seen: dict[str, set[tuple[str, str]]] = {}
    saw_model = False
    sheet_records = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record == "SHEET":
                sheet_records += 1
            if record not in {"ATOM", "HETATM"} or len(line) < 27:
                continue
            if line[12:16].strip() != "CA" or line[16].strip() not in {"", "A", "1"}:
                continue
            residue = THREE_TO_ONE.get(line[17:20].strip().upper())
            if residue is None:
                continue
            chain = line[21].strip() or "_"
            if selected_chain is not None and chain != selected_chain:
                continue
            residue_id = (line[22:26].strip(), line[26].strip())
            seen.setdefault(chain, set())
            if residue_id in seen[chain]:
                continue
            seen[chain].add(residue_id)
            residues.setdefault(chain, []).append(residue)
    for chain, sequence in residues.items():
        yield chain, "".join(sequence), {
            "format": "pdb",
            "beta_sheet_records": sheet_records,
        }


def _parse_cif_atom_loop(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read the common row-oriented _atom_site loop used by PDBx/mmCIF."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        columns: list[str] = []
        cursor = index + 1
        while cursor < len(lines) and lines[cursor].lstrip().startswith("_"):
            columns.append(lines[cursor].strip().split()[0])
            cursor += 1
        if not columns or not columns[0].startswith("_atom_site."):
            continue
        rows: list[list[str]] = []
        tokens: list[str] = []
        while cursor < len(lines):
            stripped = lines[cursor].strip()
            if not stripped or stripped.startswith("#"):
                if rows:
                    break
                cursor += 1
                continue
            if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break
            tokens.extend(shlex.split(stripped, comments=False, posix=True))
            while len(tokens) >= len(columns):
                rows.append(tokens[:len(columns)])
                tokens = tokens[len(columns):]
            cursor += 1
        return columns, rows
    return [], []


def _parse_cif(path: Path, selected_chain: str | None) -> Iterator[tuple[str, str, dict[str, Any]]]:
    columns, rows = _parse_cif_atom_loop(path)
    if not columns:
        return
    column = {name: index for index, name in enumerate(columns)}

    def value(row: list[str], *names: str, default: str = "") -> str:
        for name in names:
            if name in column and row[column[name]] not in {".", "?"}:
                return row[column[name]]
        return default

    residues: dict[str, list[str]] = {}
    seen: dict[str, set[tuple[str, str]]] = {}
    first_model: str | None = None
    for row in rows:
        if value(row, "_atom_site.group_PDB") not in {"ATOM", "HETATM"}:
            continue
        model = value(row, "_atom_site.pdbx_PDB_model_num", default="1")
        if first_model is None:
            first_model = model
        if model != first_model:
            continue
        if value(row, "_atom_site.auth_atom_id", "_atom_site.label_atom_id") != "CA":
            continue
        if value(row, "_atom_site.label_alt_id") not in {"", "A", "1"}:
            continue
        residue = THREE_TO_ONE.get(
            value(row, "_atom_site.auth_comp_id", "_atom_site.label_comp_id").upper()
        )
        if residue is None:
            continue
        chain = value(row, "_atom_site.auth_asym_id", "_atom_site.label_asym_id", default="_")
        if selected_chain is not None and chain != selected_chain:
            continue
        residue_id = (
            value(row, "_atom_site.auth_seq_id", "_atom_site.label_seq_id"),
            value(row, "_atom_site.pdbx_PDB_ins_code"),
        )
        seen.setdefault(chain, set())
        if residue_id in seen[chain]:
            continue
        seen[chain].add(residue_id)
        residues.setdefault(chain, []).append(residue)
    for chain, sequence in residues.items():
        yield chain, "".join(sequence), {"format": "mmcif", "beta_sheet_records": 0}


def _records_from_path(path: Path, selected_chain: str | None) -> Iterator[dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix in {".fa", ".faa", ".fasta"}:
        for record_name, sequence in _parse_fasta_text(path.read_text(), str(path)):
            yield analyze_sequence(sequence, record_name, source=str(path))
    elif suffix == ".pdb":
        for chain, sequence, structure in _parse_pdb(path, selected_chain):
            yield analyze_sequence(
                sequence, f"{path.stem}:{chain}", source=str(path),
                chain=chain, structure=structure,
            )
    elif suffix in {".cif", ".mmcif"}:
        for chain, sequence, structure in _parse_cif(path, selected_chain):
            yield analyze_sequence(
                sequence, f"{path.stem}:{chain}", source=str(path),
                chain=chain, structure=structure,
            )
    else:
        raise ValueError(f"unsupported file type: {path}")


def _directory_files(path: Path, recursive: bool, kind: str = "auto") -> list[Path]:
    iterator = path.rglob("*") if recursive else path.glob("*")
    suffixes = {
        "auto": SUPPORTED_SUFFIXES,
        "fasta": {".fa", ".faa", ".fasta"},
        "structure": {".pdb", ".cif", ".mmcif"},
    }[kind]
    return sorted(
        item for item in iterator
        if item.is_file() and item.suffix.lower() in suffixes
    )


def _stdin_records(text: str) -> Iterator[tuple[str, str]]:
    if text.lstrip().startswith(">"):
        yield from _parse_fasta_text(text, "stdin")
        return
    for index, line in enumerate(text.splitlines(), start=1):
        if line.strip():
            yield f"stdin:{index}", line.strip()


def collect_records(
    inputs: Iterable[str],
    selected_chain: str | None = None,
    recursive: bool = False,
    stdin: TextIO = sys.stdin,
    kind: str = "auto",
) -> list[dict[str, Any]]:
    """Load raw sequences, FASTA/structure files, directories, or stdin."""
    records: list[dict[str, Any]] = []
    for item in inputs:
        if item == "-":
            for name, sequence in _stdin_records(stdin.read()):
                records.append(analyze_sequence(sequence, name, source="stdin"))
            continue
        path = Path(item)
        if path.is_dir():
            files = _directory_files(path, recursive, kind)
            if not files:
                raise ValueError(f"no {kind} input files found in directory: {path}")
            for file_path in files:
                records.extend(_records_from_path(file_path, selected_chain))
        elif path.is_file():
            records.extend(_records_from_path(path, selected_chain))
        else:
            records.append(analyze_sequence(
                item, f"sequence_{len(records) + 1}", source="argument"
            ))
    return records


SUMMARY_HEADERS = (
    "id", "chain", "length", "hydrophobic_run", "beta_branched_run",
    "aspartimide", "cysteine", "risk",
)


def _summary_row(record: dict[str, Any]) -> dict[str, Any]:
    features = record["features"]
    return {
        "id": record["id"],
        "chain": record["chain"] or "",
        "length": features["length"],
        "hydrophobic_run": features["max_hydrophobic_run"],
        "beta_branched_run": features["max_beta_branched_run"],
        "aspartimide": features["aspartimide_motif_count"],
        "cysteine": features["cysteine_count"],
        "risk": record["risk_level"],
    }


def render_table(records: list[dict[str, Any]], output: TextIO) -> None:
    rows = [_summary_row(record) for record in records]
    widths = {header: len(header) for header in SUMMARY_HEADERS}
    for row in rows:
        for header in SUMMARY_HEADERS:
            widths[header] = max(widths[header], len(str(row[header])))
    output.write("  ".join(header.ljust(widths[header]) for header in SUMMARY_HEADERS) + "\n")
    output.write("  ".join("-" * widths[header] for header in SUMMARY_HEADERS) + "\n")
    for row in rows:
        output.write(
            "  ".join(str(row[header]).ljust(widths[header]) for header in SUMMARY_HEADERS) + "\n"
        )


def render_text(records: list[dict[str, Any]], output: TextIO) -> None:
    for record_index, record in enumerate(records):
        if record_index:
            output.write("\n")
        features = record["features"]
        output.write(f"{record['id']}  [{record['risk_level'].upper()}]\n")
        output.write(f"Sequence: {record['sequence']}\n")
        output.write(
            "Key: "
            f"length={features['length']}, "
            f"hydrophobic_run={features['max_hydrophobic_run']}, "
            f"beta_branched_run={features['max_beta_branched_run']}, "
            f"aspartimide={features['aspartimide_motif_count']}, "
            f"cysteine={features['cysteine_count']}\n"
        )
        output.write(
            f"Context: GRAVY={features['gravy']:+.3f}, "
            f"charge@7.4={features['net_charge_ph7_4']:+.2f}, "
            f"pI={features['isoelectric_point']:.2f}, "
            f"MW={features['molecular_weight_da']:.1f} Da\n"
        )
        output.write("Recommendations:\n")
        for recommendation in record["recommendations"]:
            output.write(f"- {recommendation}\n")


def render_records(records: list[dict[str, Any]], output_format: str, output: TextIO) -> None:
    if output_format == "text":
        render_text(records, output)
    elif output_format == "table":
        render_table(records, output)
    elif output_format == "json":
        json.dump(records, output, indent=2)
        output.write("\n")
    elif output_format == "jsonl":
        for record in records:
            output.write(json.dumps(record, separators=(",", ":")) + "\n")
    elif output_format == "csv":
        writer = csv.DictWriter(output, fieldnames=SUMMARY_HEADERS)
        writer.writeheader()
        writer.writerows(_summary_row(record) for record in records)
    else:
        raise ValueError(f"unsupported output format: {output_format}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract peptide SPPS features from sequences, FASTA, PDB, or mmCIF.",
        epilog="Examples: spps_check.py ACDEFG | spps_check.py designs/ --format jsonl",
    )
    parser.add_argument("inputs", nargs="*", help="Sequence, file, directory, or '-' for stdin.")
    parser.add_argument("--chain", help="Only this structure chain ID.")
    parser.add_argument("--recursive", action="store_true", help="Search directories recursively.")
    parser.add_argument(
        "--kind", choices=("auto", "fasta", "structure"), default="auto",
        help="Filter directory inputs by file kind (default: auto/all supported files).",
    )
    parser.add_argument(
        "--format", choices=("text", "table", "json", "jsonl", "csv"),
        help="Output format (default: text for one record, table for many).",
    )
    parser.add_argument(
        "--fail-on", choices=("none", "high", "medium", "any"), default="none",
        help="Exit 2 if a record reaches this risk level.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s schema {SCHEMA_VERSION}")

    # Compatibility with the original prototype CLI.
    parser.add_argument("--input_fasta_dir", help=argparse.SUPPRESS)
    parser.add_argument("--input_structure_dir", help=argparse.SUPPRESS)
    parser.add_argument("--display_table", action="store_true", help=argparse.SUPPRESS)
    return parser


def _should_fail(records: list[dict[str, Any]], fail_on: str) -> bool:
    if fail_on == "none":
        return False
    threshold = {"any": 1, "medium": 2, "high": 3}[fail_on]
    return any(SEVERITY_RANK.get(record["risk_level"], 0) >= threshold for record in records)


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not args.inputs and not args.input_fasta_dir and not args.input_structure_dir:
        parser.error("provide at least one sequence, file, directory, or '-'")

    try:
        records = collect_records(
            args.inputs, args.chain, args.recursive, kind=args.kind
        )
        if args.input_fasta_dir:
            records.extend(collect_records(
                [args.input_fasta_dir], args.chain, args.recursive, kind="fasta"
            ))
        if args.input_structure_dir:
            records.extend(collect_records(
                [args.input_structure_dir], args.chain, args.recursive, kind="structure"
            ))
        if not records:
            raise ValueError("no peptide records found (check the selected chain and input)")
        output_format = args.format or (
            "table" if args.display_table or len(records) > 1 else "text"
        )
        render_records(records, output_format, sys.stdout)
        return 2 if _should_fail(records, args.fail_on) else 0
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
