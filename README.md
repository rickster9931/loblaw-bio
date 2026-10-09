# Loblaw Bio - miraclib cell count analysis

Loads `cell-count.csv` into a SQLite database, runs the analyses for parts 2-4, and shows the results in a Streamlit dashboard.

## Running it

In a Codespace (or locally with Python 3.10+):

```
make setup       # install requirements
make pipeline    # build cell_counts.db and write results to outputs/
make dashboard   # start the dashboard on port 8501
```

Dashboard: http://localhost:8501 once `make dashboard` is running (in Codespaces, open it from the Ports tab).

`make clean` deletes the database and outputs.

## Files

- `load_data.py` - creates `cell_counts.db` and loads the csv (part 1)
- `analysis.py` - queries and stats for parts 2-4, saves results to `outputs/`
- `dashboard.py` - Streamlit app, uses the functions in `analysis.py`
- `Makefile` - setup / pipeline / dashboard targets

## Schema

```
projects    (project_id)
subjects    (subject_id, project_id, condition, age, sex, treatment, response)
samples     (sample_id, subject_id, sample_type, time_from_treatment_start)
populations (population_id, name)
cell_counts (sample_id, population_id, count)
```

In the csv, subject info like condition, treatment and response is repeated on every sample, so I moved it into a `subjects` table. Sample type and timepoint stay with the sample. Cell counts are stored one row per sample/population rather than one column per cell type, so adding a new population later just means adding rows.

For more projects and samples this should still hold up. There are indexes on the columns used for filtering, and new kinds of measurements could go in their own tables keyed by `sample_id`. If it got too big for SQLite the same schema would work in Postgres.

## Results

**Part 2** - `outputs/part2_frequency_table.csv`, one row per sample and population (52,500 rows).

**Part 3** - Melanoma patients on miraclib, PBMC samples only (993 responder samples, 975 non-responder). I used a Mann-Whitney U test for each population since the percentages aren't normally distributed, and adjusted the p-values with Benjamini-Hochberg since there are 5 tests.

| population | median diff (R - NR) | p | q |
|---|---|---|---|
| b_cell | -0.36 | 0.056 | 0.139 |
| cd8_t_cell | +0.13 | 0.639 | 0.639 |
| cd4_t_cell | +0.56 | 0.013 | 0.067 |
| nk_cell | -0.29 | 0.121 | 0.202 |
| monocyte | -0.33 | 0.163 | 0.204 |

CD4 T cells are the only population with p < 0.05, slightly higher in responders. It doesn't hold up after correcting for multiple tests (q = 0.067), and the difference is small, so I'd treat it as something to follow up on rather than a predictor.

Each subject has samples from days 0, 7 and 14, which aren't independent, so I also ran the test on per-subject averages (CD4 still p = 0.012) and on day 0 only (nothing significant). Those are in `outputs/part3_stats_*.csv` and in the dashboard. Boxplot: `outputs/part3_boxplot.png`.

**Part 4** - Melanoma PBMC samples at day 0 from miraclib patients: 656 samples.

- Samples per project: prj1 384, prj3 272
- Subjects: 331 responders, 325 non-responders
- Subjects: 344 male, 312 female

Average B cell count for melanoma males who responded, at day 0: 10206.15
