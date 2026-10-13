# analysis/

Collects the results of all jobs in `workspace/` into one table and makes standard plots.
Only reads the small text files in each job folder; no `.gsd` files and no simulation needed.

## Run

Inside Ubuntu (WSL), from the repository folder:

```bash
conda activate thiolene
python analysis/collect_results.py      # -> analysis/results/
python analysis/make_plots.py           # -> analysis/plots/
```

For another workspace (for example your own test jobs):

```bash
python analysis/collect_results.py --workspace ~/thiolene-test/workspace --out ~/thiolene-test/results
python analysis/make_plots.py --results ~/thiolene-test/results --out ~/thiolene-test/plots
```

## Output

| File | Contents |
| --- | --- |
| `results/results.xlsx` | **Read me** (every column explained), **Recipes** (mean and std over replicas), **Jobs** (one row per job) |
| `results/jobs.csv`, `results/recipes.csv` | The same two tables as CSV |
| `results/curves.csv` | Per job and saved frame: molecule conversion, C=C and SH conversion, molecule sizes |
| `plots/gel_point.png` | C=C conversion at the gel point vs tetrathiol %, with ideal Flory–Stockmayer theory |
| `plots/final_conversion.png` | Final C=C and SH conversion |
| `plots/network_structure.png` | Gel fraction, load-carrying fraction, primary loops, load-carrying strand density |
| `plots/modulus.png` | Young's modulus vs tetrathiol %, and vs load-carrying strand density |
| `plots/conversion_curves.png` | Conversion over time and molecule growth, 50 % tetrathiol |

## Things to know

- **Two conversion scales.** `trajectory_conversion.txt` and the group's gel point use *molecule conversion*
  (1 − molecules now / molecules at start). The table also gives the gel point as C=C conversion
  (`gel_point_ene_conversion`), which is what the plots use.
- An empty cell means that analysis step was never run for that job (for example, SH conversion is
  missing for the 30 stiffness-10 jobs and the 4 jobs with period 50 or 500).
- `propagation_prob` is listed as set, but the simulation always used 1.0 (known bug).
- The job document's `N_total` is wrong; `n_beads` in the table is recounted correctly.
- All values are in simulation units.
