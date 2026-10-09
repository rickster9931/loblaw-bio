"""Creates cell_counts.db and loads cell-count.csv into it.

Run with: python load_data.py
"""

import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_FILE = ROOT / "cell-count.csv"
DB_FILE = ROOT / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# subject info (condition, treatment, response, age, sex) repeats on every
# sample row in the csv, so it gets its own table. counts are stored one row
# per sample/population so adding a new cell type doesn't need a new column.
SCHEMA = """
CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition  TEXT NOT NULL,
    age        INTEGER,
    sex        TEXT CHECK (sex IN ('M', 'F')),
    treatment  TEXT NOT NULL,
    response   TEXT CHECK (response IN ('yes', 'no'))
);

CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER
);

CREATE TABLE populations (
    population_id INTEGER PRIMARY KEY,
    name          TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id     TEXT    NOT NULL REFERENCES samples(sample_id),
    population_id INTEGER NOT NULL REFERENCES populations(population_id),
    count         INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subjects_filter ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_filter ON samples(sample_type, time_from_treatment_start);
CREATE INDEX idx_samples_subject ON samples(subject_id);
"""


def to_int(value):
    return int(value) if value.strip() else None


def main():
    # start fresh every run
    if DB_FILE.exists():
        DB_FILE.unlink()

    conn = sqlite3.connect(DB_FILE)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA)

    for name in POPULATIONS:
        conn.execute("INSERT INTO populations (name) VALUES (?)", (name,))
    pop_ids = dict(conn.execute("SELECT name, population_id FROM populations"))

    with open(CSV_FILE, newline="") as f:
        for row in csv.DictReader(f):
            conn.execute("INSERT OR IGNORE INTO projects VALUES (?)", (row["project"],))
            conn.execute(
                "INSERT OR IGNORE INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    row["subject"],
                    row["project"],
                    row["condition"],
                    to_int(row["age"]),
                    row["sex"] or None,
                    row["treatment"],
                    row["response"] or None,
                ),
            )
            conn.execute(
                "INSERT INTO samples VALUES (?, ?, ?, ?)",
                (
                    row["sample"],
                    row["subject"],
                    row["sample_type"],
                    to_int(row["time_from_treatment_start"]),
                ),
            )
            for name in POPULATIONS:
                conn.execute(
                    "INSERT INTO cell_counts VALUES (?, ?, ?)",
                    (row["sample"], pop_ids[name], int(row[name])),
                )

    conn.commit()
    n_samples = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    n_subjects = conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
    conn.close()
    print(f"Loaded {n_samples} samples ({n_subjects} subjects) into {DB_FILE.name}")


if __name__ == "__main__":
    main()
