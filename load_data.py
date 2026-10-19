#!/usr/bin/env python3
"""Part 1: initialize the SQLite database and load cell-count.csv into it.

Usage:
    python load_data.py

Creates ``cell_counts.db`` in the repository root (re-created from scratch on
every run so the pipeline is idempotent).
"""

import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
PRAGMA foreign_keys = ON;

-- One row per project / study.
CREATE TABLE projects (
    project_id      TEXT PRIMARY KEY
);

-- One row per subject (patient). Attributes that are fixed for a subject
-- (demographics, indication, treatment arm, response) live here, not on samples.
CREATE TABLE subjects (
    subject_id      TEXT PRIMARY KEY,
    project_id      TEXT NOT NULL REFERENCES projects(project_id),
    condition       TEXT NOT NULL,              -- indication, e.g. melanoma / carcinoma / healthy
    age             INTEGER,
    sex             TEXT CHECK (sex IN ('M', 'F')),
    treatment       TEXT NOT NULL,              -- e.g. miraclib / phauximab / none
    response        TEXT CHECK (response IN ('yes', 'no'))  -- NULL when not applicable
);

-- One row per biological sample.
CREATE TABLE samples (
    sample_id                   TEXT PRIMARY KEY,
    subject_id                  TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                 TEXT NOT NULL,  -- e.g. PBMC / WB
    time_from_treatment_start   INTEGER NOT NULL
);

-- Lookup table of measured immune cell populations.
CREATE TABLE cell_populations (
    population_id   INTEGER PRIMARY KEY,
    name            TEXT NOT NULL UNIQUE
);

-- Long-format counts: one row per (sample, population). Adding a new
-- population later only requires a new cell_populations row, not a schema change.
CREATE TABLE cell_counts (
    sample_id       TEXT NOT NULL REFERENCES samples(sample_id),
    population_id   INTEGER NOT NULL REFERENCES cell_populations(population_id),
    count           INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subjects_filter ON subjects(condition, treatment, response);
CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_samples_filter ON samples(sample_type, time_from_treatment_start);
CREATE INDEX idx_counts_population ON cell_counts(population_id);
"""


def _clean(value):
    """Convert empty CSV fields to NULL."""
    value = value.strip()
    return value if value else None


def load(csv_path=CSV_PATH, db_path=DB_PATH):
    if db_path.exists():
        db_path.unlink()

    conn = sqlite3.connect(db_path)
    try:
        conn.executescript(SCHEMA)
        conn.executemany(
            "INSERT INTO cell_populations (population_id, name) VALUES (?, ?)",
            list(enumerate(POPULATIONS, start=1)),
        )
        pop_ids = {name: i for i, name in enumerate(POPULATIONS, start=1)}

        projects, subjects, samples, counts = set(), {}, [], []
        with open(csv_path, newline="") as fh:
            for row in csv.DictReader(fh):
                row = {k: _clean(v) for k, v in row.items()}
                projects.add(row["project"])

                subject = (
                    row["subject"], row["project"], row["condition"],
                    int(row["age"]) if row["age"] else None,
                    row["sex"], row["treatment"], row["response"],
                )
                if subjects.setdefault(row["subject"], subject) != subject:
                    raise ValueError(f"Inconsistent metadata for subject {row['subject']}")

                samples.append((
                    row["sample"], row["subject"], row["sample_type"],
                    int(row["time_from_treatment_start"]),
                ))
                counts.extend(
                    (row["sample"], pop_ids[pop], int(row[pop])) for pop in POPULATIONS
                )

        conn.executemany("INSERT INTO projects VALUES (?)", [(p,) for p in sorted(projects)])
        conn.executemany("INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)", subjects.values())
        conn.executemany("INSERT INTO samples VALUES (?, ?, ?, ?)", samples)
        conn.executemany("INSERT INTO cell_counts VALUES (?, ?, ?)", counts)
        conn.commit()
    finally:
        conn.close()

    print(
        f"Loaded {len(samples)} samples, {len(subjects)} subjects, "
        f"{len(projects)} projects, {len(counts)} cell counts into {db_path.name}"
    )


if __name__ == "__main__":
    load()
