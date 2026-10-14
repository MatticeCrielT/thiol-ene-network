# thiol-ene-network

Coarse-grained molecular dynamics of a radical thiol-ene photopolymerization, and analysis of the
network that forms (gel point, loops, load-carrying strands, modulus).

Molecules are short bead chains (a di-ene, a dithiol and a tetrathiol). HOOMD-blue moves the beads;
every `polymerize_period` timesteps the code checks which radicals touch a partner and applies three
reactions: thiol-ene propagation, chain transfer and ene homopolymerization. signac / signac-flow
manage the parameter sweep: every combination of settings is one "job" folder in `workspace/`.

**Where is everything?** See [MAP.md](MAP.md): every folder, every script, every step and every
result file.

## Setup (Windows: inside WSL / Ubuntu; HOOMD does not run on native Windows)

```bash
conda env create -f environment.yml     # once
conda activate thiolene
```

## Run

```bash
python init.py                          # make job folders from statepoint_grid in init.py
python project.py status                # what is done, what can run next
python project.py run -o equilibrate    # build and relax the mixture  -> equi.gsd
python project.py run -o polymerize     # react until conversion stops -> polymerize.gsd
python project.py run -o calculate_conversion   # then any analysis step (list in MAP.md)
```

Add `-j <job id>` to run one job, `-n 1` to run one at a time.
The existing `workspace/` jobs have 10,000 monomers (about 62,500 beads): too large for a laptop.
For a local test, set `"N_monomers": [500]` and a single replica in `init.py` (a 500-monomer run takes
about 30 minutes on a laptop CPU).

## Collect and plot results

```bash
python analysis/collect_results.py      # all jobs -> analysis/results/results.xlsx (+ csv)
python analysis/make_plots.py           # -> analysis/plots/*.png
```

See [analysis/README.md](analysis/README.md).

## Known limitations

- `environment.yml` uses HOOMD 4.9.1, which runs every step. The group's `workspace/` results were
  made with HOOMD 4.0.1, which cannot run `deform` / `deform_gel`.
- The stretching test stops at about 35-45 % elongation (a FENE bond reaches its maximum length),
  as in the group's runs. In small systems (a few thousand beads) the modulus from the first 10 %
  strain is dominated by noise; use 10,000 monomers or fit a longer part of the curve.
- Bugs found while reading the code are listed in [MAP.md](MAP.md#known-issues).

The original notes of the group are in [old/README_original](old/README_original).
