# Surface solar irradiance variability: code and processed data

This repository contains the processed data and Python code used to generate
Figs. 2–6 of the associated manuscript. The observations were collected by a
dense network of nine pyranometers over an approximately 10 km × 10 km domain
on the Inner Mongolia Plateau.

The repository starts from the processed clear-sky index used in the manuscript.
It does not include the original MATLAB data, original MATLAB scripts,
raw-data conversion utilities, or internal validation programs.

## Contents

| Path | Description |
|---|---|
| `analysis.py` | Numerical calculations for temporal variability, spatial correlation, and point-to-area differences |
| `run_all.py` | Generates Figs. 2–6 and their CSV source tables |
| `requirements.txt` | Python dependencies |
| `data/analysis_inputs.npz` | Analysis-ready clear-sky-index arrays |
| `data/clear_sky_index.csv.gz` | Portable text copy of the clear-sky-index data |
| `data/sites.csv` | Coordinates of the nine stations |
| `data/selected_days.csv` | Days selected for the three sky categories |
| `data/legacy_fits.json` | Archived fit coefficients required by Fig. 6 |

## Run

Python 3.11 or newer is recommended.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python run_all.py --mode clean --layout double --fig6-networks nine
```

The command creates `outputs/clean_double/` and writes Figs. 2–6 as 600 dpi
PNG, vector PDF, and editable SVG files. It also exports a CSV source table for
each figure and a JSON run report.

The plotting presets use 7 pt text at final size. The single-column layout is
8.5 cm wide, and the two-column layout is 17.5 cm wide. Use `--width-cm`,
`--font-pt`, and `--dpi` if the target journal specifies different values.

## Processed data

`analysis_inputs.npz` contains:

- `kt`: clear-sky index, shape `(540 minutes, 9 stations, 91 days)`;
- `site`: latitude and longitude, shape `(2, 9)`;
- `dates`: 91 ISO calendar dates;
- `cl`, `ov`, and `bc`: zero-based indices for clear-sky, overcast, and
  shallow-cumulus days.

The study subset comprises 12 clear-sky days, 8 overcast days, and 18
shallow-cumulus days. The key `bc` is retained for compatibility with the
historical analysis files. `NaN` denotes a missing station-day.

`clear_sky_index.csv.gz` contains the same clear-sky-index observations in a
portable text format. The file has one row per date and sample, followed by the
nine station values.

## Figures

| Figure | Analysis |
|---|---|
| Fig. 2 | Example daily clear-sky-index series under three sky conditions |
| Fig. 3 | Effect of temporal averaging on a shallow-cumulus day |
| Fig. 4 | Temporal variability as a function of averaging scale |
| Fig. 5 | Station-pair correlation as a function of distance/time scale |
| Fig. 6 | Difference between the nine-site area mean and the central station |

The public result uses `--mode clean`. Fig. 6 contains only the nine-site area
mean used in the revised manuscript. The right axes in Figs. 4 and 6 are
reference scales obtained by multiplying the dimensionless clear-sky-index
statistic by 300 W m−2; they are not independently calculated GHI statistics.

## Citation and license

Add the final paper citation, repository DOI, and approved code/data licenses
before making the repository public.
