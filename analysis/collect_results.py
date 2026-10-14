"""Collect the results of every job in a signac workspace into one table.

Reads the small text files that the simulation and analysis steps leave in each
job folder (no .gsd files needed) and writes:

    results/jobs.csv      one row per job: its settings + its key results
    results/recipes.csv   one row per recipe (all settings except replica_index):
                          mean, standard deviation and number of replicas
    results/curves.csv    conversion and molecule size per saved frame, per job
    results/results.xlsx  Read me + Recipes + Jobs sheets

A result that a job does not have (step never run for that job) is left empty.

Usage (inside the thiolene conda environment):
    python analysis/collect_results.py
    python analysis/collect_results.py --workspace ~/thiolene-test/workspace --out my_results
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent

# Settings from signac_statepoint.json, renamed to plain names.
SETTINGS = {
    "crosslinker_percent": "crosslinker_percent",
    "chain_side_reaction_probability": "homopolymerization_prob",
    "chain_transfer_probability": "chain_transfer_prob",
    "thiol_reaction_probability": "propagation_prob",
    "radical_percent": "radical_percent",
    "angle_constant": "angle_constant",
    "polymerize_period": "polymerize_period",
    "N_monomers": "N_monomers",
    "density": "density",
    "temperature": "temperature",
    "r_cut_reaction": "r_cut_reaction",
    "monomer_size": "monomer_size",
    "extender_size": "extender_size",
}

# Every result column with its meaning. Used for the "Read me" sheet.
COLUMNS = [
    ("job_id", "Folder name of the job in workspace/"),
    ("replica_index", "Repeat number of the same recipe (different random start)"),
    ("crosslinker_percent", "% of thiol groups that sit on the tetrathiol (rest on the dithiol)"),
    ("homopolymerization_prob", "Chance per check that C• adds to a C=C (setting chain_side_reaction_probability)"),
    ("chain_transfer_prob", "Chance per check that C• takes H from an SH"),
    ("propagation_prob", "Setting thiol_reaction_probability. NOTE: not passed to the code, the simulation always used 1.0"),
    ("radical_percent", "% of thiols turned into S• at the start"),
    ("angle_constant", "Bond-angle stiffness (higher = stiffer chains)"),
    ("polymerize_period", "Timesteps between reaction checks"),
    ("N_monomers", "Number of ene monomers (each has 2 C=C)"),
    ("density", "Beads per unit volume at the start"),
    ("temperature", "kT in simulation units"),
    ("r_cut_reaction", "Distance within which beads can react (bead diameters)"),
    ("monomer_size", "Spacer setting (0 = simple 4-bead ene monomer)"),
    ("extender_size", "Spacer setting (0 = simple 2-bead dithiol)"),
    ("n_tetrathiols", "Number of tetrathiol molecules"),
    ("n_dithiols", "Number of dithiol molecules"),
    ("n_beads", "Total beads in the box (simple molecules only)"),
    ("total_timesteps", "Last timestep of the polymerization (polymerize.txt)"),
    ("final_ene_conversion", "Fraction of C=C groups reacted at the end (job document)"),
    ("final_thiol_conversion", "Fraction of SH groups reacted at the end (trajectory_thiol_ene_conversion.txt)"),
    ("ene_used_by_homopolymerization", "Share of reacted C=C that reacted with another C=C instead of a thiol = (ene - thiol conversion) / ene conversion"),
    ("final_molecule_conversion", "1 - molecules now / molecules at start (the scale the group used for the gel point)"),
    ("volume_shrinkage_percent", "Box volume decrease from start to end of the polymerization, % (NPT at zero pressure)"),
    ("gel_point_molecule_conversion", "Gel point: molecule conversion where the reduced weight-average size peaks (gelation2.py)"),
    ("gel_point_ene_conversion", "Gel point expressed as C=C conversion (same frame)"),
    ("gel_point_thiol_conversion", "Gel point expressed as SH conversion (same frame)"),
    ("gel_point_frame", "Saved frame of polymerize.gsd at the gel point"),
    ("gel_fraction", "Beads in the largest molecule / all beads, at the end"),
    ("n_strands", "Linear strands in the gel between branch points (strand_sizes.txt)"),
    ("mean_strand_length", "Number-average strand length in beads"),
    ("n_dangling_strands", "Dangling ends hanging off the gel (loop_counts.txt)"),
    ("dangling_beads", "Beads in dangling ends"),
    ("dangling_bead_fraction", "Dangling-end beads / beads in the gel"),
    ("n_primary_loops", "Loops of 1 strand: a strand that starts and ends on the same junction"),
    ("n_secondary_loops", "Loops of 2 strands between the same 2 junctions"),
    ("n_tertiary_loops", "Loops of 3 strands"),
    ("n_quaternary_loops", "Loops of 4 strands"),
    ("n_primary_loops_tetrathiol", "Primary loops closing on a tetrathiol (primary_loop_types.txt)"),
    ("n_primary_loops_ene", "Primary loops closing on an ene carbon (primary_loop_types.txt)"),
    ("elastically_effective_fraction", "Beads that carry load / beads in the gel (contract_bonds_analysis/overall.txt)"),
    ("n_junctions", "Beads in the gel with 3 or more bonds (branch points)"),
    ("n_effective_junctions", "Junctions with 3+ load-carrying strands (Scanlan)"),
    ("n_bridging_junctions", "Junctions with exactly 2 load-carrying strands (act as a chain link)"),
    ("n_ineffective_junctions", "Junctions with 0-1 load-carrying strands"),
    ("avg_effective_functionality", "Average number of load-carrying strands per effective junction"),
    ("n_effective_junctions_tetrathiol", "Effective junctions that are a tetrathiol centre"),
    ("n_effective_junctions_carbon", "Effective junctions that are a carbon made by homopolymerization"),
    ("n_effective_strands", "Load-carrying strands (effective_strand_hist.txt)"),
    ("mean_effective_strand_length", "Number-average length of load-carrying strands, beads"),
    ("effective_strand_density", "Load-carrying strands per unit volume (job document crosslinking_density)"),
    ("youngs_modulus", "Young's modulus of the isolated gel, simulation units (youngs_modulus_gel.txt)"),
]


# ---------------------------------------------------------------------------
# Small readers. Each returns None (or empty) when the file is missing.
# ---------------------------------------------------------------------------
def read_number(path):
    try:
        return float(path.read_text().split()[0])
    except (FileNotFoundError, IndexError, ValueError):
        return None


def read_table(path, skip_header=1, delimiter=None):
    """Numeric rows of a text file as a 2D array, skipping the header lines."""
    if not path.exists():
        return None
    rows = []
    for line in path.read_text().splitlines()[skip_header:]:
        parts = line.replace(",", " ").split() if delimiter == "," else line.split()
        try:
            rows.append([float(p) for p in parts])
        except ValueError:
            continue
    return np.array(rows) if rows else None


def read_header_row(path):
    """Files with one header line and one value line, tab-separated -> dict."""
    if not path.exists():
        return {}
    lines = path.read_text().splitlines()
    return {k: float(v) for k, v in zip(lines[0].split(), lines[1].split())}


def read_loop_counts(path):
    """loop_counts.txt holds five (size count) blocks, each after a '... data' line."""
    blocks = {}
    if not path.exists():
        return blocks
    current = None
    for line in path.read_text().splitlines():
        if line.strip().endswith("data"):
            current = line.split()[0].split("_")[0]  # dangling / primary / ...
            blocks[current] = []
        elif current and line.strip():
            size, count = line.split()
            blocks[current].append((float(size), float(count)))
    return blocks


def read_polymerize_log(path):
    """Last timestep and the first/last box volume from polymerize.txt."""
    if not path.exists():
        return None, None
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if parts and parts[0].isdigit():
            rows.append(parts)
    if not rows:
        return None, None
    return int(rows[-1][0]), float(rows[-1][6])


def count_effective_junction_types(path):
    if not path.exists():
        return None, None
    table = pd.read_csv(path, sep="\t")
    effective = table[table["scanlan_classification"] == "effective"]
    return int((effective["type"] == "Tetra-S").sum()), int((effective["type"] == "Carbon").sum())


# ---------------------------------------------------------------------------
# One job -> one row (+ its curves)
# ---------------------------------------------------------------------------
def collect_job(folder):
    sp = json.loads((folder / "signac_statepoint.json").read_text())
    doc = json.loads((folder / "signac_job_document.json").read_text())
    row = {"job_id": folder.name, "replica_index": sp.get("replica_index")}
    for key, name in SETTINGS.items():
        row[name] = sp.get(key)

    # molecule counts (the job document's N_total is wrong, so recount here)
    row["n_tetrathiols"] = doc.get("N_crosslinkers")
    row["n_dithiols"] = doc.get("N_extenders")
    simple = sp.get("monomer_size", 0) == 0 and sp.get("extender_size", 0) == 0
    if simple and row["n_tetrathiols"] is not None:
        row["n_beads"] = 5 * row["n_tetrathiols"] + 2 * row["n_dithiols"] + 4 * sp["N_monomers"]

    # --- reaction ----------------------------------------------------------
    row["total_timesteps"], final_volume = read_polymerize_log(folder / "polymerize.txt")
    if final_volume and row.get("n_beads"):
        start_volume = row["n_beads"] / sp["density"]
        row["volume_shrinkage_percent"] = 100 * (1 - final_volume / start_volume)

    row["final_ene_conversion"] = doc.get("reacted_monomers")
    molecules = read_table(folder / "trajectory_conversion.txt")       # frame, conv, largest, avg, Mw
    thiol_ene = read_table(folder / "trajectory_thiol_ene_conversion.txt")  # frame, thiol, ene
    if molecules is not None:
        row["final_molecule_conversion"] = molecules[-1, 1]
        if row.get("n_beads"):
            row["gel_fraction"] = molecules[-1, 2] / row["n_beads"]
    if thiol_ene is not None:
        row["final_thiol_conversion"] = thiol_ene[-1, 1]
        if row["final_ene_conversion"]:
            row["ene_used_by_homopolymerization"] = (
                row["final_ene_conversion"] - row["final_thiol_conversion"]) / row["final_ene_conversion"]

    # --- gel point ---------------------------------------------------------
    gel = doc.get("gelation_conversion_2_2")
    if gel:
        row["gel_point_molecule_conversion"], row["gel_point_frame"] = gel[0], int(gel[1])
        if thiol_ene is not None:
            at_gel = thiol_ene[thiol_ene[:, 0] == int(gel[1])]
            if len(at_gel):
                row["gel_point_thiol_conversion"], row["gel_point_ene_conversion"] = at_gel[0, 1], at_gel[0, 2]

    # --- network structure at the end ---------------------------------------
    strands = read_table(folder / "strand_sizes.txt")
    if strands is not None:
        row["n_strands"] = strands[:, 1].sum()
        row["mean_strand_length"] = (strands[:, 0] * strands[:, 1]).sum() / strands[:, 1].sum()

    loops = read_loop_counts(folder / "loop_counts.txt")
    if loops:
        dangling = loops.get("dangling", [])
        row["n_dangling_strands"] = sum(c for _, c in dangling)
        row["dangling_beads"] = sum(s * c for s, c in dangling)
        if row.get("gel_fraction") and row.get("n_beads"):
            row["dangling_bead_fraction"] = row["dangling_beads"] / (row["gel_fraction"] * row["n_beads"])
        for kind in ["primary", "secondary", "tertiary", "quaternary"]:
            row[f"n_{kind}_loops"] = sum(c for _, c in loops.get(kind, []))

    primary_types = read_header_row(folder / "primary_loop_types.txt")
    row["n_primary_loops_tetrathiol"] = primary_types.get("N_Tetrathiol_primary_loops")
    row["n_primary_loops_ene"] = primary_types.get("N_ene_primary_loops")

    # --- elastic network (after contract_bonds) ------------------------------
    cb = folder / "contract_bonds_analysis"
    overall = read_header_row(cb / "overall.txt")
    if overall:
        eff, ineff = overall["n_effective_particles"], overall["n_ineffective_particles"]
        row["elastically_effective_fraction"] = eff / (eff + ineff)
    scanlan = read_header_row(cb / "scanlan_case_analysis.txt")
    if scanlan:
        row["n_junctions"] = scanlan["n_crosslinkers"]
        row["n_effective_junctions"] = scanlan["n_effective_crosslinkers"]
        row["n_bridging_junctions"] = scanlan["n_bridging_crosslinkers"]
        row["n_ineffective_junctions"] = scanlan["n_ineffective_crosslinkers"]
        row["avg_effective_functionality"] = scanlan["avg_effective_Functionality"]
    row["n_effective_junctions_tetrathiol"], row["n_effective_junctions_carbon"] = \
        count_effective_junction_types(cb / "crosslink_properties.txt")
    effective = read_table(cb / "effective_strand_hist.txt", skip_header=0, delimiter=",")  # count, size
    if effective is not None and effective[:, 0].sum() > 0:
        row["n_effective_strands"] = effective[:, 0].sum()
        row["mean_effective_strand_length"] = (effective[:, 0] * effective[:, 1]).sum() / effective[:, 0].sum()
    row["effective_strand_density"] = doc.get("crosslinking_density")

    row["youngs_modulus"] = read_number(folder / "youngs_modulus_gel.txt")

    # --- curves --------------------------------------------------------------
    curves = None
    if molecules is not None:
        curves = pd.DataFrame(molecules[:, :5], columns=[
            "frame", "molecule_conversion", "largest_molecule", "average_molecule", "weight_average_molecule"])
        curves["frame"] = curves["frame"].astype(int)
        if thiol_ene is not None:
            te = pd.DataFrame(thiol_ene[:, :3], columns=["frame", "thiol_conversion", "ene_conversion"])
            te["frame"] = te["frame"].astype(int)
            curves = curves.merge(te, on="frame", how="left")
        # a frame is saved every 10 reaction checks
        curves.insert(1, "timestep", curves["frame"] * 10 * sp["polymerize_period"])
        curves.insert(0, "job_id", folder.name)
    return row, curves


def summarize_recipes(jobs):
    """Mean, standard deviation and replica count per recipe."""
    settings = [c for c in SETTINGS.values() if c in jobs.columns]
    results = [c for c in jobs.columns if c not in settings + ["job_id", "replica_index"]]
    grouped = jobs.groupby(settings, dropna=False)[results]
    summary = grouped.mean().add_suffix("_mean").join(grouped.std().add_suffix("_std"))
    summary.insert(0, "n_replicas", jobs.groupby(settings, dropna=False).size())
    # put mean and std of each result next to each other
    ordered = ["n_replicas"] + [f"{c}_{s}" for c in results for s in ("mean", "std")]
    return summary[ordered].reset_index()


def write_excel(path, jobs, recipes):
    readme = pd.DataFrame(COLUMNS, columns=["column", "meaning"])
    readme = pd.concat([pd.DataFrame([
        {"column": "Recipes sheet", "meaning": "One row per recipe: _mean and _std over its replicas, n_replicas = number of repeats"},
        {"column": "Jobs sheet", "meaning": "One row per simulation; empty cell = that analysis was not run for that job"},
        {"column": "Units", "meaning": "Simulation (reduced) units: lengths in bead diameters, energy in kT"},
        {"column": "", "meaning": ""},
    ]), readme])
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, frame in [("Read me", readme), ("Recipes", recipes), ("Jobs", jobs)]:
            frame.to_excel(xl, sheet_name=name, index=False)
            sheet = xl.sheets[name]
            sheet.freeze_panes = "B2"
            for col in sheet.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col[:50])
                sheet.column_dimensions[col[0].column_letter].width = min(max(width + 2, 10), 60 if name == "Read me" else 28)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workspace", type=Path, default=REPO / "workspace", help="signac workspace folder")
    parser.add_argument("--out", type=Path, default=REPO / "analysis" / "results", help="output folder")
    args = parser.parse_args()

    folders = sorted(p for p in args.workspace.iterdir() if (p / "signac_statepoint.json").exists())
    rows, curves = [], []
    for folder in folders:
        row, curve = collect_job(folder)
        rows.append(row)
        if curve is not None:
            curves.append(curve)

    order = [c for c, _ in COLUMNS]
    jobs = pd.DataFrame(rows)
    jobs = jobs[[c for c in order if c in jobs.columns]]
    jobs = jobs.sort_values(["angle_constant", "polymerize_period", "homopolymerization_prob",
                             "crosslinker_percent", "replica_index"]).reset_index(drop=True)
    recipes = summarize_recipes(jobs)

    args.out.mkdir(parents=True, exist_ok=True)
    jobs.to_csv(args.out / "jobs.csv", index=False)
    recipes.to_csv(args.out / "recipes.csv", index=False)
    if curves:
        pd.concat(curves).to_csv(args.out / "curves.csv", index=False)
    write_excel(args.out / "results.xlsx", jobs, recipes)

    print(f"{len(jobs)} jobs, {len(recipes)} recipes -> {args.out}")
    filled = jobs.notna().sum()
    missing = filled[filled < len(jobs)]
    if len(missing):
        print("Columns not available for every job (jobs with a value):")
        for col, n in missing.items():
            print(f"  {col}: {n}")


if __name__ == "__main__":
    main()
