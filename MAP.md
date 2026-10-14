# Repository map

Status of 2026-10-14 (branch `mattice/bugfixes`). What every folder and file is, how the steps connect,
and what every result file means.

## 1. Top level

```
thiol-ene-network/
├── README.md              start here: what it is, setup, how to run
├── MAP.md                 this file
├── environment.yml        conda environment (Linux / WSL)
├── init.py                YOU EDIT: the settings grid -> one job folder per combination
├── project.py             the list of steps (operations) and when each is done (labels)
├── scripts/               simulation and analysis code used by project.py
├── analysis/              collect all results into one table + standard plots (added 2026-10)
├── workspace/             the group's 124 finished jobs (text results only, no .gsd)
├── percolation-analyzer/  compiled tool that writes percolation_data.txt
├── signac_project_document.json   signac bookkeeping (operation status), do not edit
├── .signac/               marks this folder as a signac project, do not edit
└── old/                   earlier phases, tests and one-off helpers (not needed to run)
```

## 2. The steps and what they write

Run each with `python project.py run -o <operation>`. A step can run when its "needs" condition
holds and is skipped when its "done when" condition holds. `.gsd` files are large and are **not
in git**; the `workspace/` folders in git only hold the text results.

**Simulation** (code in `scripts/simulate.py` + `scripts/system.py` unless noted)

| Operation | Needs | Done when / writes | What it does |
|---|---|---|---|
| `equilibrate` | - | `equi.gsd` | Build molecules, remove overlaps, 5,000 steps at kT = 1 |
| `polymerize` | `equi.gsd` | `polymerize.gsd`, `polymerize.txt`, job doc `reacted_monomers` | Radicals on, reaction checks + motion until conversion stops |
| `contract_bonds` | reacted | `contract_bonds.gsd` | Shrink every bond to ~0 so loops and dangling ends collapse; used to find load-carrying strands |
| `deform_gel` | reacted | `deform_gel.gsd`, `deform_gel.log`, `pressures_gel.gsd` | Stretch the isolated gel (`scripts/deform_gel.py`). Needs HOOMD >= 4.6 |
| `deform` | reacted | `deform.gsd`, `deform.log` | Older stretch of the whole box (`scripts/deform.py`). Needs HOOMD >= 4.6 |
| `diffuse` | `equi.gsd` | `diffuse.gsd`, `diffuse.txt` | Motion without reactions; only for 50 %, homopolymerization 0, replica 0 |
| `percolation_contract_bonds` | `percolation_frame` in job doc | `percolation_contract_bonds.gsd` | `contract_bonds` on the frame stored as `percolation_frame` |

"reacted" = `calculate_conversion` has run after `polymerize` and C=C conversion has not risen over
the last `stop_after_rounds` reaction checks (molecule conversion if the C=C file is missing).

**How `polymerize` ends.** It stops when conversion has not risen for `stop_after_rounds` reaction
checks (job document; `init.py` sets 50 for ~200 radicals and proportionally more for fewer radicals,
e.g. 1000 for 500 monomers), or after 10,000 checks per run (1 million timesteps at period 100). In
the second case run `polymerize` again: it continues from `polymerize.gsd`. Repeat
`polymerize` + `calculate_conversion` until `project.py status` shows the job as reacted.
Analyses use the last saved frame and are not refreshed when a polymerization is continued later:
delete their output files first.

**Analysis** (code in `scripts/network_properties.py` unless noted)

| Operation | Needs | Writes | Content |
|---|---|---|---|
| `calculate_conversion` | `polymerize.gsd` | `trajectory_conversion.txt`, `trajectory_thiol_ene_conversion.txt` | Per frame: molecule conversion, largest / average / weight-average molecule size; SH and C=C conversion (`scripts/conversion.py`) |
| `gelation_2_analysis` | reacted | job doc `gelation_conversion_2_2` = [molecule conversion, frame] | Gel point = peak of the weight-average size without the largest molecule (`scripts/gelation2.py`) |
| `strand_analysis` | reacted | `strand_sizes.txt` (+ `file.txt`) | Strand lengths in the gel between branch points; `file.txt` = gel fraction (bug: meant to be `gel_fraction.txt`) |
| `defect_analysis` | reacted | `loop_counts.txt` | Dangling ends and primary to quaternary loops, as size-count lists |
| `primary_loop_types_analysis` | reacted | `primary_loop_types.txt` | Primary loops on a tetrathiol vs. on an ene carbon |
| `xlink_rdf_analysis` | reacted | `thiol_thiol_rdf.txt`, `ene_ene_rdf.txt`, `thiol_ene_rdf.txt` | Radial distribution of crosslink points (how evenly they are spread) |
| `conv95_defect_analysis`, `conv95_strand_analysis` | reacted | `loop_counts_conversion95.txt`, `conv95_strand_sizes.txt` | Same as above at the first frame where molecule conversion passes 0.95 |
| `contract_bonds_analysis` | `contract_bonds.gsd` | `contract_bonds_analysis/` 5 files (below) | Load-carrying (elastically effective) beads, strands and junctions |
| `calculate_cld` | contract analysis | job doc `crosslinking_density` | Load-carrying strands per volume. Its "done" check always says no, so it always shows as runnable |
| `vv_analysis` | contract analysis | `contract_bonds_analysis/voronoi_volumes_{all,thiol,ene}.txt` | Voronoi volume around crosslinks (local packing) |
| `deformation_gel_analysis` | `deform_gel` | `stress_strain_gel.txt` | Columns: true strain, stress, uniaxial strain (lambda^2 - 1/lambda), tensile stress (`scripts/plot_deforms.py`) |
| `modulus_gel_analysis` | stress-strain | `youngs_modulus_gel.txt` | Young's modulus = 3 x slope of tensile stress vs. uniaxial strain, first 10 % strain |
| `diffusion_analysis` | `diffuse` | `diffusion_length_at_delta_timesteps_<period>.txt` | How far molecules move in one reaction period (`scripts/diffusion.py`) |
| `get_percolation_frame` | `percolation_data.txt` | job doc `percolation_frame` | Meant: frame where the gel first spans the box. Actually records the last such frame (see known issues) (`scripts/conversion.py`) |
| `percolation_contract_bonds_analysis`, `calculate_percolation_cld` | percolation contract | `percolation_contract_bonds_analysis/`, job doc `perc_crosslinking_density` | The contract analysis at the percolation frame |

`percolation_data.txt` is not made by `project.py`: it comes from
`percolation-analyzer/bin/sample_graph -i polymerize.gsd -o percolation_data.txt` (see
`old/helper_scripts/percolate.py`). The tool's source code is not in the repository; the compiled
Linux program runs in WSL.

## 3. `scripts/` (all used by project.py)

| File | Role |
|---|---|
| `system.py` | The chemistry: bead / bond types, `molecule_counts()`, building the molecules (spacer length via `monomer_size` / `extender_size`), initiation and the three reactions |
| `simulate.py` | The physics: force field, `equilibrate`, `polymerize` (reaction loop), `contract_bonds`, `diffuse` |
| `customAction/polymerize.py` | Alternative GPU reaction code for `polymerization_method = custom_action_*`; not used by the current jobs (`cpu_local_snapshot`) |
| `conversion.py` | Conversion and molecule sizes per frame; 95 % and percolation frames |
| `extract.py` | Helper (`connected_components`) used by `conversion.py` |
| `gelation2.py` | Gel point |
| `network_properties.py` | Strands, loops, dangling ends, RDFs, contract-bond analysis, Scanlan junction classes, Voronoi volumes |
| `deform_gel.py`, `deform.py` | Stretching tests |
| `plot_deforms.py` | Stress-strain curve and Young's modulus from the stretching test |
| `diffusion.py` | Diffusion analysis |

## 4. Result files in each job folder (`workspace/<job id>/`)

| File | Made by | Content |
|---|---|---|
| `signac_statepoint.json` | `init.py` | The settings of this job (do not edit) |
| `signac_job_document.json` | several | Counts, `reacted_monomers` (final C=C conversion), gel point, crosslink densities, timings. `N_total` and `L` in it are wrong (miscounted beads) |
| `polymerize.txt` | `polymerize` | Log every 10 checks: timestep, speed, temperature, pressure, energies, box volume |
| `trajectory_conversion.txt` | `calculate_conversion` | frame, molecule conversion, largest, average and weight-average molecule size (beads) |
| `trajectory_thiol_ene_conversion.txt` | `calculate_conversion` | frame, SH conversion, C=C conversion |
| `strand_sizes.txt`, `loop_counts.txt`, `primary_loop_types.txt` | analyses above | Size-count lists |
| `*_rdf.txt` | `xlink_rdf_analysis` | r, g(r) |
| `contract_bonds_analysis/overall.txt` | contract analysis | Load-carrying vs. non-load-carrying beads in the gel |
| `contract_bonds_analysis/effective_strand_hist.txt`, `ineffective_strand_hist.txt` | contract analysis | count, strand length |
| `contract_bonds_analysis/scanlan_case_analysis.txt` | contract analysis | Junctions: total, effective (3+ load-carrying strands), bridging (2), ineffective; average functionality |
| `contract_bonds_analysis/crosslink_properties.txt` | contract analysis | Per junction: bead id, type (Tetra-S or Carbon), effective functionality, class |
| `stress_strain_gel.txt`, `youngs_modulus_gel.txt` | stretching analysis | See section 2 |
| `percolation_data.txt` | percolation-analyzer | frame, percolation dimension (0-3), largest cluster, number of molecules |
| `gel_fraction.txt`, `crosslink_density_1.txt`, `cluster_sizes.txt`, `dangle_sizes.txt`, `molecule_size_histogram.json`, `trajectory_molecule_size.txt`, `stress_strain.txt`, `youngs_modulus.txt` | **older code** (no current step makes them) | Kept for reference; `gel_fraction.txt` = beads in largest molecule / all beads |

Two conversion scales are used: *molecule conversion* = 1 - molecules now / molecules at start (the
group's gel point and `trajectory_conversion.txt`), and *C=C / SH conversion* = share of groups
reacted. `analysis/` gives the gel point on both.

## 5. Results available now

| Where | What |
|---|---|
| `workspace/` (124 jobs) | Main study: 10,000 monomers; tetrathiol 0-100 % in steps of 12.5; homopolymerization 0 and 0.15; stiffness 100 (all) and 10 (12.5 / 50 / 100 %); 5 replicas each; period 100 (+ 50 and 500 for one recipe). Not every analysis was run for every job; `analysis/collect_results.py` lists which |
| `analysis/results/results.xlsx` | All 124 jobs and 28 recipes in one table, every column explained on the "Read me" sheet |
| `analysis/plots/` | Gel point, final conversion, network structure, modulus, conversion curves |
| `~/thiolene-test` (in Ubuntu, not in git) | 500-monomer test run of 2026-10-12, including its `.gsd` movie |
| `old/archive_workspace/` | Earlier phase (79 jobs, 2024-25): different physics (density 0.9, kT 0.9, stiffness 5), mostly 500 monomers |

## 6. `old/`

| Folder / file | What it was |
|---|---|
| `README_original` | The group's original notes (requirements and physical parameters describe the earlier phase) |
| `project_old.py` | Earlier version of `project.py` (was the file `test`) |
| `archive_workspace/`, `archive_workspace_test/` | Earlier simulation phase and its tests |
| `performance_testing/` | Speed tests of CPU vs GPU reaction code, with copies of `system.py` / `simulate.py` |
| `plots/` | The group's figures from 2024-25, mostly of the earlier phase |
| `analyze/` | 2024 scripts to draw the network as a graph (netgraph) |
| `scripts/` | Stand-alone plot / helper scripts that `project.py` does not use, and `deform_test/` images |
| `helper_scripts/` | One-off root scripts: run all equilibrations, percolation runs for 4 jobs, progress overview, renaming files, restarting failed jobs. They expect to be run from the repository root and contain hard-coded job ids or paths |
| `job_progress/` | Output of `helper_scripts/progress.py` (2024) |
| `config/freud.ini` | Settings of an unrelated HTTP tool, not of the freud library |

## Known issues

Still open:

| Issue | Where | Effect |
|---|---|---|
| `deform` / `deform_gel` use `hoomd.variant.box` | `scripts/deform*.py` | Needs HOOMD >= 4.6; `environment.yml` now uses 4.9.1 (tested) |
| Stretching stops with "Error in bond calculation" at strain 0.3-0.4 (FENE bond at its limit) | `scripts/deform_gel.py` | Same in the group's runs; the stress-strain file is still written and used |
| An interrupted `deform_gel` leaves files that count as "done" | `project.py` label `deformed_gel` | Delete `deform_gel.gsd`, `deform_gel.log`, `pressures_gel.gsd` and run again |
| `calculate_cld` done-check always false | `project.py` | Always listed as runnable |
| 34 of the group's jobs have no `trajectory_thiol_ene_conversion.txt` | `workspace/` | Cannot be recreated without their `.gsd` files (fixed for new jobs) |
| 8 of the 60 thiol-ene-only jobs (stiffness 100) end at C=C conversion 0.924-0.9255 | `workspace/` | Points to a stop at 92.5 % in an earlier code version (`scripts/simulate.py` still has it for the GPU path); their final values may be slightly low |
| With the stricter `reacted` label (C=C conversion), 9 group jobs count as not finished | `workspace/` | They rose < 0.7 % over their last 100 frames; results usable, but their `.gsd` files are not here to continue |
| Job `03f1d195...` (12.5 %) has a wrong `N_crosslinkers` (1250) in its job document | `workspace/` | The simulation used 625 (checked from its first frame); `analysis/` recounts from the recipe |
| Homopolymerization may add to the other C=C of the same monomer (a code comment says this should be excluded) | `scripts/system.py` `chain_growth` | Rare: that C=C is usually further away than the reaction distance |

Fixed on 2026-10-14 (branch `mattice/bugfixes`). Results in `workspace/` were made **before** these
fixes; only the percolation frame changes their meaning:

| Issue | Fix |
|---|---|
| `thiol_reaction_probability` never reached the reaction code (always 1.0) | Passed on in `scripts/simulate.py`. All existing jobs used 1, so they are unaffected |
| Spaced molecules could never switch on; their builder was broken | `monomer_size` / `extender_size` = number of spacer beads (type D) in the di-ene / dithiol. Size 0 gives bit-identical systems to the old builder (tested) |
| Bead counts in `init.py` differed from `system.py` (`N_total`, `L` wrong in job documents) | Both use `molecule_counts()` in `scripts/system.py`. Old job documents still hold the wrong values |
| `get_percolation_frame` stored the last spanning frame instead of the first | Stops at the first. The 20 existing "percolation" analyses describe the final network, not the gel point |
| Gel fraction written to `file.txt` | Writes `gel_fraction.txt` again |
| `percolation_analyzed` label returned nothing | Returns its value |
| Neighbour filter used `|` instead of `&` | Fixed (no change in behaviour) |
| `calculate_conversion` never created `trajectory_thiol_ene_conversion.txt` (only updated it) | Creates it for new jobs |
| Load-carrying analysis broke for dithiols with spacers (spacer bonds labelled `-SH` set the cutoff: 1 % instead of 87 % load-carrying) | Only `-SH` / C=C bonds on unreacted beads set the cutoff; spacer bonds labelled C-C single. Unchanged for simple molecules (tested) |
| The `reacted` label used molecule conversion over 4 frames, so unfinished small jobs counted as done | Uses C=C conversion over `stop_after_rounds` (same 4 frames for the default 50) |
| `polymerize` always stopped after 50 checks without a rise; small systems stopped far too early (a 500-monomer job at 58 %) | Window `stop_after_rounds` in the job document, set by `init.py` from the number of radicals; default stays 50 |
