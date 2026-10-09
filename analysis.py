"""Analysis for parts 2-4. Reads from cell_counts.db (run load_data.py first).

Running this file saves all the results to outputs/. dashboard.py imports
the same functions so both show the same numbers.
"""

import sqlite3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # no GUI needed, we only save the figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
DB_FILE = ROOT / "cell_counts.db"
OUT_DIR = ROOT / "outputs"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
ALPHA = 0.05
COLORS = {"Responder": "#2a78d6", "Non-responder": "#eb6834"}


def connect(check_same_thread=True):
    if not DB_FILE.exists():
        raise FileNotFoundError(f"{DB_FILE.name} not found, run load_data.py first.")
    return sqlite3.connect(DB_FILE, check_same_thread=check_same_thread)


# ---- Part 2 ----

FREQ_SQL = """
SELECT
    cc.sample_id AS sample,
    SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
    p.name AS population,
    cc.count AS count,
    100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
FROM cell_counts cc
JOIN populations p ON p.population_id = cc.population_id
ORDER BY cc.sample_id, p.population_id
"""


def frequency_table(conn):
    return pd.read_sql_query(FREQ_SQL, conn)


# ---- Part 3 ----

# melanoma + miraclib + PBMC, only samples where we know the response
RESPONDERS_SQL = f"""
WITH freq AS ({FREQ_SQL})
SELECT f.sample, f.population, f.percentage,
       s.subject_id, s.time_from_treatment_start, sub.response
FROM freq f
JOIN samples s ON s.sample_id = f.sample
JOIN subjects sub ON sub.subject_id = s.subject_id
WHERE sub.condition = 'melanoma'
  AND sub.treatment = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND sub.response IN ('yes', 'no')
"""


def responder_frequencies(conn):
    df = pd.read_sql_query(RESPONDERS_SQL, conn)
    df["group"] = df["response"].map({"yes": "Responder", "no": "Non-responder"})
    return df


def compare_groups(df):
    # Mann-Whitney since the percentages aren't normally distributed,
    # then BH correction because we're testing 5 populations at once
    rows = []
    for pop in POPULATIONS:
        d = df[df["population"] == pop]
        yes = d.loc[d["response"] == "yes", "percentage"]
        no = d.loc[d["response"] == "no", "percentage"]
        u, p = stats.mannwhitneyu(yes, no, alternative="two-sided")
        rows.append({
            "population": pop,
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responders": np.median(yes),
            "median_non_responders": np.median(no),
            "median_diff": np.median(yes) - np.median(no),
            "U": u,
            # effect size, > 0 means higher in responders
            "rank_biserial": 2 * u / (len(yes) * len(no)) - 1,
            "p_value": p,
        })
    res = pd.DataFrame(rows)
    res["q_value"] = stats.false_discovery_control(res["p_value"], method="bh")
    res["significant_nominal"] = res["p_value"] < ALPHA
    res["significant"] = res["q_value"] < ALPHA
    return res


def responder_statistics(df):
    # Each subject has up to 3 samples (day 0/7/14) which aren't independent,
    # so also check per-subject averages and baseline-only as sanity checks.
    per_subject = df.groupby(["subject_id", "population", "response"], as_index=False)["percentage"].mean()
    baseline = df[df["time_from_treatment_start"] == 0]
    return {
        "all_samples": compare_groups(df),
        "subject_mean": compare_groups(per_subject),
        "baseline": compare_groups(baseline),
    }


def plot_boxplots(df, path):
    groups = ["Responder", "Non-responder"]
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(16, 4.5))
    for ax, pop in zip(axes, POPULATIONS):
        data = [df.loc[(df["population"] == pop) & (df["group"] == g), "percentage"] for g in groups]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, medianprops={"color": "black"})
        for box, g in zip(bp["boxes"], groups):
            box.set_facecolor(COLORS[g])
            box.set_alpha(0.75)
        ax.set_xticks([1, 2], groups)
        ax.set_title(pop)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Relative frequency (%)")
    fig.suptitle("Melanoma, miraclib, PBMC: responders vs non-responders")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---- Part 4 ----

BASELINE_SQL = """
SELECT s.sample_id, s.subject_id, sub.project_id, sub.response, sub.sex
FROM samples s
JOIN subjects sub ON sub.subject_id = s.subject_id
WHERE sub.condition = 'melanoma'
  AND sub.treatment = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND s.time_from_treatment_start = 0
ORDER BY s.sample_id
"""


def baseline_subset(conn):
    return {
        "samples": pd.read_sql_query(BASELINE_SQL, conn),
        "samples_per_project": pd.read_sql_query(
            f"SELECT project_id AS project, COUNT(*) AS n_samples "
            f"FROM ({BASELINE_SQL}) GROUP BY project_id ORDER BY project_id", conn),
        "subjects_by_response": pd.read_sql_query(
            f"SELECT response, COUNT(DISTINCT subject_id) AS n_subjects "
            f"FROM ({BASELINE_SQL}) GROUP BY response ORDER BY response", conn),
        "subjects_by_sex": pd.read_sql_query(
            f"SELECT sex, COUNT(DISTINCT subject_id) AS n_subjects "
            f"FROM ({BASELINE_SQL}) GROUP BY sex ORDER BY sex", conn),
    }


def avg_b_cells_melanoma_male_responders_t0(conn):
    # any sample type / treatment here, not just the baseline subset above
    avg = conn.execute("""
        SELECT AVG(cc.count)
        FROM cell_counts cc
        JOIN populations p ON p.population_id = cc.population_id
        JOIN samples s ON s.sample_id = cc.sample_id
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE p.name = 'b_cell'
          AND sub.condition = 'melanoma'
          AND sub.sex = 'M'
          AND sub.response = 'yes'
          AND s.time_from_treatment_start = 0
    """).fetchone()[0]
    return round(avg, 2)


def main():
    OUT_DIR.mkdir(exist_ok=True)
    conn = connect()

    # part 2
    freq = frequency_table(conn)
    freq.to_csv(OUT_DIR / "part2_frequency_table.csv", index=False)
    print(f"Part 2: frequency table ({len(freq)} rows)")
    print(freq.head(10).to_string(index=False))

    # part 3
    resp = responder_frequencies(conn)
    plot_boxplots(resp, OUT_DIR / "part3_boxplot.png")
    results = responder_statistics(resp)
    print("\nPart 3: responders vs non-responders")
    cols = ["population", "n_responders", "n_non_responders", "median_diff",
            "rank_biserial", "p_value", "q_value", "significant_nominal", "significant"]
    for name, table in results.items():
        table.to_csv(OUT_DIR / f"part3_stats_{name}.csv", index=False)
        print(f"\n[{name}]")
        print(table[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    main_res = results["all_samples"]
    sig = main_res.loc[main_res["significant"], "population"].tolist()
    nominal = main_res.loc[main_res["significant_nominal"], "population"].tolist()
    print(f"\nSignificant after BH correction (q < {ALPHA}): {', '.join(sig) or 'none'}")
    print(f"Significant before correction (p < {ALPHA}): {', '.join(nominal) or 'none'}")

    # part 4
    subset = baseline_subset(conn)
    for name, table in subset.items():
        table.to_csv(OUT_DIR / f"part4_{name}.csv", index=False)
    print(f"\nPart 4: {len(subset['samples'])} baseline melanoma/miraclib/PBMC samples")
    for name in ["samples_per_project", "subjects_by_response", "subjects_by_sex"]:
        print(f"\n{name}:\n{subset[name].to_string(index=False)}")

    avg_b = avg_b_cells_melanoma_male_responders_t0(conn)
    (OUT_DIR / "part4_avg_b_cell.txt").write_text(f"{avg_b:.2f}\n")
    print(f"\nAverage B cells (melanoma, male, responders, day 0): {avg_b:.2f}")

    conn.close()


if __name__ == "__main__":
    main()
