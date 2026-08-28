# MWT Tap-Habituation Data Processing Pipeline

**Operations Manual — Steps 1 through 5**

This document describes the end-to-end pipeline that converts raw Multi-Worm Tracker (MWT) recordings into the summarized, statistically annotated tables that back the MWT Dashboard. It consolidates the guidance that currently lives in individual notebook cells into one reference.

---

## Table of Contents

1. [Pipeline at a Glance](#1-pipeline-at-a-glance)
2. [Prerequisites](#2-prerequisites)
3. [Data and Directory Conventions](#3-data-and-directory-conventions)
4. [Step 1 — Choreography Extraction](#4-step-1--choreography-extraction)
5. [Step 2 — Tap Speed (Per-Worm) Analysis](#5-step-2--tap-speed-per-worm-analysis)
6. [Step 3 — Tap Response Analysis](#6-step-3--tap-response-analysis)
7. [Step 4 — Baseline and Post-Stimulus Analysis](#7-step-4--baseline-and-post-stimulus-analysis)
8. [Step 5 — Statistics and Database Upload](#8-step-5--statistics-and-database-upload)
9. [Parameter Reference](#9-parameter-reference)
10. [File Reference](#10-file-reference)
11. [Database Table Reference](#11-database-table-reference)
12. [Troubleshooting](#12-troubleshooting)
13. [Full-Run Checklist](#13-full-run-checklist)

---

## 1. Pipeline at a Glance

```
Raw MWT recordings (.zip)
         │
         ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 1 — Choreography_script.ipynb                          │
│ Runs Chore.jar twice over every .zip                        │
│   Pass A (plate-level)  → .trv, .prv, .txt, plate .dat      │
│   Pass B (-N all)       → per-worm .dat  (*.00065.dat)      │
└─────────────────────────────────────────────────────────────┘
         │                    │                    │
    per-worm .dat        plate .dat             .trv
         │                    │                    │
         ▼                    ▼                    ▼
┌──────────────────┐  ┌──────────────────┐  ┌──────────────────┐
│ STEP 2           │  │ STEP 4           │  │ STEP 3           │
│ Tap speed        │  │ Baseline + PSA   │  │ Tap response     │
│                  │  │                  │  │                  │
│ *_tapspeed_      │  │ *_baseline_      │  │ needs Step 2's   │
│   byplate.csv ───┼──┼──────────────────┼─▶│ byplate.csv      │
│                  │  │   output.csv     │  │                  │
│                  │  │ *_post_          │  │ *_tap_output.csv │
│                  │  │   stimulus.csv   │  │                  │
└──────────────────┘  └──────────────────┘  └──────────────────┘
                              │                    │
                              └────────┬───────────┘
                                       ▼
                    ┌─────────────────────────────────────┐
                    │ STEP 5 — data_processing_for_db     │
                    │ Merge → derive metrics → t-tests →  │
                    │ FDR → MSD → PostgreSQL upload       │
                    └─────────────────────────────────────┘
                                       │
                                       ▼
                              MWT Dashboard (Streamlit)
```

**Mandatory execution order:** Step 1 → Step 2 → Step 3 → Step 4 → Step 5.

Steps 2 and 4 are independent of each other and can be run in either order, but **Step 3 will fail without Step 2's output**, because it merges the `1s_Speed` column from `{Screen}_tapspeed_byplate.csv`. Step 5 requires all three upstream CSVs to be present in the same folder.

---

## 2. Prerequisites

### Software

| Requirement | Notes |
|---|---|
| Java runtime | Needed for `Chore.jar` in Step 1 |
| `Chore.jar` | Must sit in the same directory the notebook is launched from |
| Python 3.11+ | Kernel used by the notebooks |
| Jupyter Lab / Notebook | Widgets must be enabled for `FileChooser` to render |

### Python packages

```
pandas  numpy  scipy  statsmodels  seaborn  matplotlib
ipywidgets  ipyfilechooser  tqdm
sqlalchemy  psycopg  requests
```

Step 5 additionally imports a local module, `backend_config`, which must be importable from the notebook's working directory.

### Hardware

Step 1 launches the JVM with `-Xms13g`. The machine needs at least 16 GB of RAM, and realistically more if anything else is running. Lower this flag if you are working on a smaller machine.

Steps 2 and 4 load every `.dat` row into memory before aggregating. On large screens this is the memory bottleneck of the pipeline. Both notebooks contain float64 → float32 downcasting; Step 2 applies it automatically, Step 4 has it commented out and available if needed.

### Database credentials

Step 5 reads connection settings via `load_config()` from a `database.ini` file. Populate `user`, `password`, `host`, `port`, and `database` before running the upload cells. The notebook exits if `user` or `password` is blank.

---

## 3. Data and Directory Conventions

**This section is the single most common source of pipeline failures. Read it before running anything.**

### Required folder hierarchy

Every notebook derives metadata by splitting the absolute file path. That means the directory tree is not cosmetic — it is the metadata schema.

```
<any parent folders>/
└── <Screen>/                          ← path[-4]  e.g. PD_Screen
    └── <Gene_Allele>/                 ← path[-3]  e.g. lrrk2_tm1898, or N2
        └── <YYYYMMDD_HHMMSS>/         ← path[-2]  e.g. 20100422_112007
            ├── <plate>.zip
            ├── <plate>.trv
            ├── <plate>.dat
            └── <plate>.00065.dat ...
```

Derived fields:

| Field | Derivation | Example |
|---|---|---|
| `Screen` | 4th path element from the end | `PD_Screen` |
| `dataset` | 3rd path element from the end | `lrrk2_tm1898` |
| `Gene` | `dataset` split on the first `_` | `lrrk2` |
| `Allele` | remainder after the first `_`; `N2` if absent | `tm1898` |
| `Date` | experiment folder name before the first `_` | `20100422` |
| `Plate_id` | `<experiment folder>_<last underscore token of filename>` | `20100422_112007_E0418aa` |
| `worm_id` | numeric token before `.dat` (Step 2 only) | `00065` |

`Plate_id` is the real join key across all steps. The `plate` column is only an ordinal within a strain and is dropped or ignored at every merge.

### Strain matching is substring matching

`ProcessData()` selects files with `strain in filepath` — a plain substring test against the **entire path**. Consequences:

- If any folder anywhere in the selected path contains a strain name, every file underneath it will be attributed to that strain. A parent folder named `N2_experiments` will make the whole tree look like N2.
- Strain names that are substrings of other strain names will over-match.

**Fix:** rename intermediate folders so they do not contain strain tokens. Use dates, or lower-case names, for any folder above the `<Gene_Allele>` level.

### macOS AppleDouble files

Files beginning with `._` are skipped explicitly in all steps. They are harmless but will inflate file counts in the "number of files found" printouts.

### Control strains

Every statistical comparison is made against N2. Two conventions exist:

- **Most screens:** the control genotype folder is `N2`.
- **`Neuron_Genes_Screen`:** two controls, `N2_N2` and `N2_XJ1`. Steps 2 and 4 reorder the strain dictionary to put both first; Step 5 branches on the screen name to build the correct control mask.

If a screen has no `N2` folder, Steps 2 and 4 raise a `ValueError` at the `genotypes.pop(genotypes.index("N2"))` line.

---

## 4. Step 1 — Choreography Extraction

**Notebook:** `Step1_Choreography_script.ipynb`
**Input:** `.zip` archives of raw MWT recordings
**Output:** `.trv`, `.prv`, `.txt`, and `.dat` files, written beside each `.zip`

### Purpose

Wraps Rex Kerr's `Chore.jar` and runs it over every `.zip` archive found under the selected folder. Two separate passes are required because they produce different products.

### Procedure

1. **Launch the notebook from the directory containing `Chore.jar`.** Cell 2 records `repo_directory = os.getcwd()` and builds the jar path from it. If you started Jupyter elsewhere, the subprocess call will fail silently per file.
2. Run the import cell.
3. Run the `FileChooser` cell and select the data folder. The default starting directory is `/Users`; edit `starting_directory` if your data lives on an external volume.
4. Confirm `folder_path` printed correctly.
5. Run the filelist cell. It walks the tree and collects every `.zip`. Check the printed count against what you expect.
6. **Run Pass A** (the cell under *For Analyzing Tap Habituation Statistics*).
7. Run the `.trv` check cell. It lists any `.zip` whose extracted directory contains no `.trv`. Investigate anything reported.
8. **Run Pass B** (the cell under *To analyze response speed (individual worms)*).
9. Run the per-worm `.dat` check cell.

### The two passes

Both passes share the same core flags:

```
-p 0.027        pixel size (mm)
-s 0.1          minimum time a worm must be tracked
-t 20           minimum duration for a valid object
-M 2            minimum travel distance (body lengths)
--shadowless    exclude objects in shadowed regions
-S              use a segmented/spine-based analysis
-o nNss*b12xyMmeSakcrP    output column specification
--plugin Reoutline::exp
--plugin Respine
```

**Pass A — plate-level, reversal measurement.** Adds:

```
--plugin MeasureReversal::tap::dt=1::collect=0.5::postfix=trv
--plugin MeasureReversal::puff::dt=3::collect=0.5::postfix=prv
--plugin MeasureReversal::postfix=txt
```

Produces one `.trv` per plate (tap-evoked reversals, consumed by Step 3), one `.prv` (puff-evoked), one `.txt`, and plate-aggregate `.dat` files (consumed by Step 4).

**Pass B — per-worm.** Adds `-N all` and omits the MeasureReversal plugins. Produces one `.dat` per tracked worm, named `<plate>.<5-digit worm id>.dat` (consumed by Step 2).

### Re-running failures

Both check cells populate `missed_filelist`. To retry only the failures, edit the processing cell to iterate over `missed_filelist` instead of `filelist` — the alternate loop line is already present, commented out, at the top of each processing cell.

The `.dat` check distinguishes the two passes by filename: it only counts a `.dat` as a Pass B product if the token before the extension is exactly five digits.

### Runtime

This is by far the longest step. Expect hours for a full screen. There is no resume-in-place; the check cells exist so that an interrupted run can be restarted against the remainder.

---

## 5. Step 2 — Tap Speed (Per-Worm) Analysis

**Notebook:** `Step2_TapSpeed_Screen_Data_Analysis.ipynb`
**Input:** per-worm `.dat` files from Step 1 Pass B (`*.<digits>.dat`)
**Output:**
- `{Screen}_tapspeed_processed_data.csv` — full retained per-frame data
- `{Screen}_tapspeed_byplate.csv` — plate × tap means, **required by Step 3**

### Purpose

Measures how fast individual worms reverse in the first second after each tap. This produces the `1s_Speed` metric, which is not available from the plate-level `.trv` data.

### Procedure

| Cell group | Action |
|---|---|
| 1 | Run imports |
| 2 | Select folder via `FileChooser`, then select the screen from the dropdown and run the cell that assigns `Screen` |
| 3 | Set `number_of_taps`, `ISI`, `first_tap` — verify the printed tolerance windows look right |
| 4 | Build filelist (regex `\.\d+\.dat$`) — check the count and the example path |
| 5 | Define processing functions (no input) |
| 6.1 | Build `StrainNames` dictionary — confirm the genotype count and that N2 is first |
| 6.2 | Process all strains (long-running, shows a progress bar) |
| 7 | Concatenate, label, aggregate, and export |

### The retention-window filter

`filter_retention_windows()` is the analytical core of this step and has no equivalent elsewhere in the pipeline. For each worm trace it retains only the frames that constitute a genuine tap-evoked reversal:

- **Window opens** when `Tap` transitions 0 → 1 while `Bias` is ≥ 0 (worm was not already reversing).
- **Window closes** at the earlier of: 1 second after `Bias` first reaches −1, or the frame before `Bias` transitions from negative back to non-negative.
- **Window is discarded** if `Bias` never reaches −1 within 1 second of opening — the worm did not respond.
- **Window is discarded** if `Bias` hits −1 on the same frame as the tap transition — treated as a spontaneous reversal that happened to coincide with the tap.

Worms with no retained windows are skipped, but the experiment counter still increments so numbering stays consistent.

### Aggregation

After retention filtering, rows are kept only where `Bias == -1` (actively reversing), then aggregated twice:

1. Per (tap, worm, plate): `Time` = min, `Speed` = mean.
2. Per (Screen, dataset, Gene, Allele, Date, Plate_id, plate, tap): mean of the above, with `Speed` renamed to `1s_Speed`.

The second table is what Step 3 consumes.

### Notes

- Tap tolerance windows here are wide — effectively the full ISI (599.9–609.9 s for tap 1) — because per-frame data must be attributed to whichever tap interval it falls in, not to a narrow event window.
- `base['taps'] = base['taps'].ffill()` forward-fills tap labels before export. Rows before the first tap carry tap 0.
- Float columns are downcast to float32 on read to control memory.

---

## 6. Step 3 — Tap Response Analysis

**Notebook:** `Step3_Tap_Screen_Response_Data_Analysis.ipynb`
**Input:** `.trv` files from Step 1 Pass A, **plus** `{Screen}_tapspeed_byplate.csv` from Step 2
**Output:** `{Screen}_tap_output.csv`

### Purpose

Extracts the classic tap-habituation response measures — reversal probability, duration, and distance — from the plate-level `.trv` files, then attaches the per-worm `1s_Speed` from Step 2.

### Procedure

Same shape as Step 2: imports → folder and screen selection → parameters → filelist → function definitions → per-strain processing → grouping and export.

### `.trv` column mapping

`.trv` files are space-delimited and headerless. Step 3 renames the columns it uses:

| Index | Name | Meaning |
|---|---|---|
| 0 | `time` | tap time |
| 2 | `rev_before` | already reversing at stimulus |
| 3 | `no_rev` | did not reverse |
| 4 | `stim_rev` | reversed in response |
| 7–15 | `dist`, `dist_std`, `dist_stderr`, `dist_0th`…`dist_100th` | reversal distance and distribution |
| 18–26 | `dura`, `dura_std`, `dura_stderr`, `dura_0th`…`dura_100th` | reversal duration and distribution |

Derived:

```
prob  = stim_rev / (no_rev + stim_rev)
speed = dist / dura
```

### Tap assignment

Tolerances here are tight — `first_tap ± 2 s`, stepped by ISI, plus a final window of (1188, 1191) for the recovery tap. This works because `.trv` rows are discrete events, not continuous frames.

The printed "total taps" per strain is worth watching. A strain with fewer rows than `plates × 31` is missing taps somewhere.

### The Step 2 merge

```python
merge_columns = ['Screen', 'dataset', 'Gene', 'Allele', 'Date', 'Plate_id', 'taps']
tapspeed_byplate = pd.read_csv(f'{Screen}_tapspeed_byplate.csv').drop(columns='plate')
TotalConcatenated = TotalConcatenated.merge(..., how='left', validate='many_to_one')
```

Three failure modes here:

- **`FileNotFoundError`** — Step 2 was not run, or was run with a different `Screen` selected, or its output is in a different folder.
- **`MergeError` from `validate='many_to_one'`** — the byplate file contains duplicate key combinations. Usually means Step 2 was run twice over overlapping data.
- **Silent `NaN` in `1s_Speed`** — key mismatch, most often a `Date` dtype difference or a plate present in the `.trv` set but absent from the per-worm `.dat` set. A commented diagnostic cell near the end lists alleles with missing `1s_Speed`; uncomment it if the column looks sparse.

### Note on the folder-naming warning

The `**Important**` banner in this notebook's Step 2 markdown cell is the substring-matching issue described in [Section 3](#3-data-and-directory-conventions). It applies to Steps 2, 3, and 4 equally.

---

## 7. Step 4 — Baseline and Post-Stimulus Analysis

**Notebook:** `Step4_Tap_Baseline_Data_Analyses.ipynb`
**Input:** plate-level `.dat` files from Step 1 Pass A
**Output:**
- `{Screen}_baseline_output.csv`
- `{Screen}_post_stimulus.csv`

### Purpose

Produces two distinct datasets from the same source files:

- **Baseline** — locomotion and morphology in the 490–590 s window, before any stimulus. This is the unstimulated phenotype.
- **Post-stimulus arousal (PSA)** — the same measures in a window 7–9.5 s after each tap, capturing the arousal state after the reversal itself has ended.

### File selection

```python
if name.endswith('.dat') and "_" in name.split(".")[-2]:
```

This selects **plate-level** `.dat` files (whose penultimate dot-token is the plate name, containing underscores) and excludes the per-worm files from Step 1 Pass B (whose penultimate token is a numeric worm id). This is the mirror image of Step 2's regex.

### `.dat` column mapping

| Index | Name | | Index | Name |
|---|---|---|---|---|
| 0 | `Time` | | 10 | `Morphwidth` |
| 1 | `n` | | 11 | `Midline` |
| 2 | `Number` | | 12 | `Area` |
| 3 | `Speed` | | 13 | `Angular Speed` |
| 4 | `Interval Speed` | | 14 | `Aspect Ratio` |
| 5 | `Bias` | | 15 | `Kink` |
| 6 | `Tap` | | 16 | `Curve` |
| 7 | `Puff` | | 17 | `Crab` |
| 8 | `x` | | 18 | `Pathlength` |
| 9 | `y` | | | |

### Parameter placement (important)

Unlike the other notebooks, `number_of_taps`, `ISI`, and `first_tap` are **not** set in Step 3 of this notebook — they are set in a cell inside section **6.2**, immediately above the processing loop. Section 3 here only defines `bins`. Set the tap parameters in 6.2, and note that the cell still prints `"done step 3"`.

The tolerances built there are `psa_tolerances`: `first_tap + 7.0` to `first_tap + 9.5`, stepped by ISI, with a final window of (1197.5, 1199). These deliberately sit *after* the reversal has finished, which is what makes the measure post-stimulus arousal rather than the response itself.

### Baseline construction

```python
Baseline_data = base.drop(columns=["Tap","Puff","x","y","Experiment","taps","plate"]).dropna()
Baseline_data = Baseline_data[(Baseline_data.Time <= 590.0) & (Baseline_data.Time >= 490.0)]
```

Note this is per-frame data, not aggregated. `{Screen}_baseline_output.csv` is the largest file the pipeline produces; Step 5 aggregates it to plate level.

### PSA construction

Rows after 599 s are grouped by (Experiment, Strain, Screen, Date, Plate_id, Gene, Allele, dataset, taps) and averaged, with `Time` taking the minimum rather than the mean. Because `taps` was assigned using `psa_tolerances`, only frames in the 7–9.5 s post-tap windows carry a tap number and survive the grouping.

An alternative implementation that builds windows explicitly from each `Tap == 1` event is present but commented out. The tolerance-based approach is the one in use.

### Known warnings

`assign_taps()` uses `df['taps'][0] = 0`, which triggers a pandas chained-assignment warning. It is suppressed by the `warnings.filterwarnings` cell and is harmless under current pandas versions, but it is a candidate for cleanup (`df.loc[df.index[0], 'taps'] = 0`, as Step 2 already does).

---

## 8. Step 5 — Statistics and Database Upload

**Notebook:** `Step5_data_processing_for_db.ipynb`
**Input:** `*baseline_output.csv`, `*tap_output.csv`, `*post_stimulus.csv` — all globbed from the selected folder
**Output:** eight CSVs plus eight PostgreSQL tables

This is the longest notebook and the only one that writes to shared infrastructure. Read [Section 8.6](#86-database-upload--read-before-running) before running the upload cells.

### 8.1 Load and merge

The three CSVs are located by glob pattern, so the folder must contain exactly one of each. PSA columns are prefixed `PSA ` on load to distinguish them from the tap-response columns of the same name.

`tap_output` and `psa_output` are merged (outer) on `dataset, Gene, Allele, Strain, Date, Plate_id, Screen, taps` to form `tap_psa_output`. Several PSA columns are dropped at this point as not analytically used: `Morphwidth`, `Midline`, `Area`, `Angular Speed`.

### 8.2 Derived habituation metrics

The tap data is sliced into three reference points and combined:

| Slice | Definition | Prefix |
|---|---|---|
| First tap | `taps == 1` | `init_` |
| Final taps | mean of `taps` 28–30 | `final_` |
| Recovery tap | `taps == 31` | `recov_` |

From these, per plate:

| Metric family | Formula |
|---|---|
| `habit_*` | `init − final` |
| `recovery_*` | `(recov − init) / init × 100` |
| `memory_retention_*` | `recov − final` |
| `sensitization_*` | `max − init` (PSA measures only) |
| `mean_*` | mean across all taps |
| `max_*` | max across all taps (PSA measures only) |

Applied across: `dura`, `prob`, `speed`, `1s_speed`, and the PSA measures `psa_speed`, `psa_bias`, `psa_aspect_ratio`, `psa_kink`, `psa_curve`, `psa_crab`.

**Screen-dependent NaN handling.** For `Neuron_Genes_Screen` and `G-Proteins_Screen`, missing values are dropped on a specific column subset. For all other screens, missing values are filled with 0. This branch exists because some screens lack a recovery tap; be aware that a filled 0 is not distinguishable downstream from a measured 0.

### 8.3 Mean sample distance (MSD)

`get_combined_MSD()` computes, for each phenotype and each gene (or allele):

- mean, count, and standard error across plates
- 95% confidence bounds via `scipy.stats.t.interval`
- all three values expressed **relative to N2**, by subtracting the N2 mean

Column names are expanded to human-readable form at this stage (`habit_dura` → `Habituation of Response Duration`, and so on) — this rename map is the canonical source for dashboard labels. The resulting MultiIndex is flattened with `-` joins, giving columns like `Initial Response Duration-mean`, `Initial Response Duration-ci95_hi`.

Run at both **gene** level (`by='Gene'`) and **allele** level (`by='dataset'`).

### 8.4 T-statistics

`do_TTest()` performs Welch's t-test (`equal_var=False`) of each gene or allele against N2, with one important refinement: **the comparison is date-matched.** For each test group, only N2 plates run on the same dates as that group are used as controls. This controls for day-to-day variation in worm batches and tracker conditions.

`pair_pvals()` then applies Benjamini–Hochberg FDR correction **row-wise** (across all metrics within a gene) at `alpha = 0.1`, and stores each cell as a `(t_statistic, fdr_corrected_p)` tuple.

`transform_tap_tstat_heatmap()` produces the plotting-ready version: any cell whose corrected p-value is ≥ 0.05 is zeroed out, so the dashboard heatmap shows only significant effects. Both raw (tuple) and transformed (numeric) versions are uploaded.

Run at both gene and allele level, for both baseline and tap data, giving four t-stat tables that are merged into two.

### 8.5 CSV exports

Before upload, columns are reordered to put identifiers first, and eight files are written to the selected folder:

```
{Screen}_baseline_output.csv              (overwritten, now reordered)
{Screen}_tap_psa_output.csv
{Screen}_final_tstat.csv
{Screen}_transformed_final_tstat.csv
{Screen}_final_tstat_allele.csv
{Screen}_transformed_final_tstat_allele.csv
{Screen}_combined_MSD.csv
{Screen}_allele_combined_MSD.csv
```

An *Alternative Method* cell reads these eight files back in. Use it to re-run the upload without recomputing the statistics — useful when an upload fails partway.

### 8.6 Database upload — read before running

**The upload cell currently uses `if_exists='replace'` for every table.** This drops and recreates each table, discarding all data from every other screen. The safer append-with-conflict-skip path is present in each line but commented out:

```python
# tap_psa_output.to_sql('tap_response_data', engine, if_exists='append',
#                       index=False, method=postgres_skip_on_duplicate)
tap_psa_output.to_sql('tap_response_data', engine, if_exists='replace',
                      index=False, method=None)
```

Before running this cell, decide which mode you want:

- **Rebuilding the whole database from scratch, one screen only** → `replace` is correct.
- **Adding or updating a single screen alongside existing data** → switch to the commented `append` lines. `postgres_skip_on_duplicate` uses `ON CONFLICT DO NOTHING` against the table's primary keys, so existing rows are preserved and duplicates are skipped.

Note that `append` mode *skips* conflicting rows rather than updating them. To genuinely update a screen already in the database, delete that screen's rows first, then append.

A final cell at the bottom of the notebook uploads a single table in isolation — use it to repair one table without touching the rest.

---

## 9. Parameter Reference

### Experiment constants

Set identically in Steps 2, 3, and 4 (in Step 4, inside section 6.2):

| Parameter | Default | Meaning |
|---|---|---|
| `number_of_taps` | 30 | Number of taps in the habituation train |
| `ISI` | 10 | Inter-stimulus interval, seconds |
| `first_tap` | 600 | Time of the first tap, seconds |
| `bins` | `linspace(0, 1200, 1201)` | 1-second bins across the recording (Steps 2, 4) |

Taps are numbered 1…31, where tap 31 is the recovery tap delivered after the ISI-long rest.

### Tap tolerance windows

The three notebooks build different windows from the same constants, because they measure different things. For tap *k* (1-indexed):

| Notebook | Window | Tap 1 | Recovery tap (31) | Rationale |
|---|---|---|---|---|
| Step 2 | `first_tap − 0.1 + 10(k−1)` to `first_tap + 9.9 + 10(k−1)` | 599.9 – 609.9 | 1190.9 – 1201 | Full ISI; every frame must be attributed to an interval |
| Step 3 | `first_tap − 2 + 10(k−1)` to `first_tap + 2 + 10(k−1)` | 598 – 602 | 1188 – 1191 | Tight window around a discrete event |
| Step 4 | `first_tap + 7.0 + 10(k−1)` to `first_tap + 9.5 + 10(k−1)` | 607 – 609.5 | 1197.5 – 1199 | Post-reversal arousal window |

If your protocol differs from 30 taps / 10 s ISI / 600 s first tap, **open a `.trv` file and read the actual tap times** before adjusting these, then verify the printed tolerance list in each notebook.

### Screen selector

All four analysis notebooks present the same dropdown. The selected value is stamped into the `Screen` column and into every output filename, so **it must be identical across Steps 2–5 for a given dataset.**

```
PD_Screen              ASD_Screen              G-Proteins_Screen
Glia_Genes_Screen      Neuron_Genes_Screen     ASD_WGS_Screen
PD_GWAS_Locus22_Screen PD_GWAS_Locus71_Screen  Miscellaneous
```

---

## 10. File Reference

### Produced by Step 1

| Pattern | Pass | Level | Consumed by |
|---|---|---|---|
| `<plate>.trv` | A | plate | Step 3 |
| `<plate>.prv` | A | plate | — (puff data, unused) |
| `<plate>.txt` | A | plate | — |
| `<plate>.dat` | A | plate | Step 4 |
| `<plate>.<worm>.dat` | B | worm | Step 2 |

### Produced by Steps 2–5

| File | Produced by | Consumed by |
|---|---|---|
| `{Screen}_tapspeed_processed_data.csv` | Step 2 | — (archive / QC) |
| `{Screen}_tapspeed_byplate.csv` | Step 2 | **Step 3** |
| `{Screen}_tap_output.csv` | Step 3 | **Step 5** |
| `{Screen}_baseline_output.csv` | Step 4 | **Step 5** (and rewritten by it) |
| `{Screen}_post_stimulus.csv` | Step 4 | **Step 5** |
| `{Screen}_tap_psa_output.csv` | Step 5 | database |
| `{Screen}_final_tstat[_allele].csv` | Step 5 | database |
| `{Screen}_transformed_final_tstat[_allele].csv` | Step 5 | database |
| `{Screen}_combined_MSD.csv` | Step 5 | database |
| `{Screen}_allele_combined_MSD.csv` | Step 5 | database |

Keep all of these in one folder per screen. Step 5's glob patterns assume it.

---

## 11. Database Table Reference

| Table | Source dataframe | Primary keys |
|---|---|---|
| `tap_response_data` | `tap_psa_output` | Date, Plate_id, Screen, taps, dataset, Gene, Allele |
| `tap_baseline_data` | `baseline_output` | Time, Plate_id, Date, Screen, dataset, Gene, Allele |
| `tstat_gene_data_raw` | `final_tstat` | Gene, Screen |
| `tstat_gene_data` | `transformed_final_tstat` | Gene, Screen |
| `tstat_allele_data_raw` | `final_tstat_allele` | dataset, Screen |
| `tstat_allele_data` | `transformed_final_tstat_allele` | dataset, Screen |
| `gene_MSD` | `combined_MSD` | Gene, Screen |
| `allele_MSD` | `allele_combined_MSD` | dataset, Screen |
| `Gene_Allele_WormBaseID` | — | WBGene, WBAllele — maintained separately, not written by this pipeline |

**Deprecated tables**, no longer written: `psa_summarized_data`, `allele_profile_data`, `gene_profile_data`. The code paths that produced them are commented out in Step 5.

### Note on normalization and pivoting

Two transformations that used to live in Step 5 have been moved into the dashboard: z-score normalization of the t-stat matrix, and melting into long format for profile plots. The commented blocks remain in `merge_Tstats()` as a record. Do not re-enable them without checking whether the dashboard still expects to do the work itself.

---

## 12. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Step 1 produces nothing, no error | Notebook not launched from the `Chore.jar` directory | Restart Jupyter from the repo root |
| Step 1: `OutOfMemoryError` / JVM won't start | `-Xms13g` exceeds available RAM | Lower to `-Xms4g` or as the machine allows |
| `FileNotFoundError: No target .dat files...` (Step 2) | Step 1 Pass B not run | Run the `-N all` cell in Step 1 |
| `FileNotFoundError: No .trv files found` (Step 3) | Step 1 Pass A not run | Run the MeasureReversal cell in Step 1 |
| `AssertionError: X is not a good identifier` | A strain in `StrainNames` matches zero files | Usually a folder named for a strain with no data in it; check the genotype folder list |
| Every strain reports the same file count | Substring collision — a parent folder contains a strain name | Rename the parent folder (see Section 3) |
| `ValueError` at `genotypes.index("N2")` | No N2 control folder in the tree | Add the control data, or adjust the control block for that screen |
| Step 3: `FileNotFoundError` on `tapspeed_byplate.csv` | Step 2 not run, or run with a different `Screen` | Re-run Step 2 with the matching screen selection |
| Step 3: `1s_Speed` all NaN after merge | Merge key mismatch, usually `Date` dtype or missing per-worm data for those plates | Uncomment the diagnostic cell listing alleles with missing `1s_Speed` |
| Step 3: `MergeError` on `validate='many_to_one'` | Duplicate keys in the byplate file | Re-run Step 2 on a clean folder |
| "This strain has N total taps" is lower than expected | Missing taps in the `.trv` files | Check the raw recording; the plate may be short |
| Memory exhaustion in Step 2 or 4 | Whole-screen per-frame data in RAM | Enable the float32 downcast cell in Step 4; process in sub-batches |
| Step 5: `IndexError` on a glob | A required upstream CSV is missing, or more than one matches | Ensure exactly one of each pattern in the folder |
| Step 5: "Please set your user and password" | `database.ini` incomplete | Populate the config file |
| Other screens vanished from the dashboard | Upload cell ran in `replace` mode | Switch to the commented `append` lines and re-upload the affected screens |
| Everything zero in the dashboard heatmap | Expected — `transform_tap_tstat_heatmap` zeroes non-significant cells (p ≥ 0.05) | Check the `_raw` tables for the underlying values |

### General notebook hygiene

- Run cells top to bottom. Several steps depend on state set in earlier cells (`folder_path`, `Screen`, `filelist`, `StrainNames`).
- Every step prints a `done step N` confirmation. A missing confirmation means the cell errored.
- Widget state is not the same as variable state: selecting a folder in the `FileChooser` does nothing until you run the following cell that assigns `folder_path`. Same for the screen dropdown.
- `os.chdir(folder_path)` is called in Steps 2, 3, and 4. Some outputs are written to the working directory rather than an absolute path, so they land in the data folder as a side effect of that call.

---

## 13. Full-Run Checklist

Processing one screen end to end:

- [ ] Verify the directory hierarchy: `<Screen>/<Gene_Allele>/<YYYYMMDD_HHMMSS>/`
- [ ] Verify no parent folder name contains a strain token
- [ ] Verify an N2 control folder exists (or `N2_N2` + `N2_XJ1` for `Neuron_Genes_Screen`)
- [ ] Confirm tap parameters against an actual `.trv` file

**Step 1**
- [ ] Launch Jupyter from the `Chore.jar` directory
- [ ] Run Pass A; confirm zero missed `.trv` directories
- [ ] Run Pass B; confirm zero missed per-worm `.dat` directories

**Step 2**
- [ ] Select folder and screen
- [ ] Verify tolerance printout
- [ ] Confirm `.dat` file count is plausible
- [ ] Confirm genotype count and that N2 sorts first
- [ ] Confirm `{Screen}_tapspeed_byplate.csv` was written

**Step 3**
- [ ] Select the **same** folder and screen
- [ ] Confirm `.trv` count
- [ ] Watch per-strain tap totals for shortfalls
- [ ] Confirm the merge succeeded and `1s_Speed` is populated
- [ ] Confirm `{Screen}_tap_output.csv` was written

**Step 4**
- [ ] Select the same folder and screen
- [ ] Set tap parameters in section **6.2**, not section 3
- [ ] Confirm both `{Screen}_baseline_output.csv` and `{Screen}_post_stimulus.csv` were written

**Step 5**
- [ ] Confirm all three input CSVs are present, exactly one match each
- [ ] Run through the statistics sections; spot-check `final_tstat.head()`
- [ ] Confirm the eight CSVs were written
- [ ] **Decide `replace` vs `append` before running the upload cell**
- [ ] Verify row counts in the database after upload
- [ ] Spot-check the dashboard against the source CSVs

---

*Maintained alongside the MWT Dashboard operations documentation. When notebook logic changes, update the affected section here in the same commit.*
