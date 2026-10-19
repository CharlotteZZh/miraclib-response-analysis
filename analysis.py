"""Analysis functions for Parts 2-4. All data is read from the SQLite database
created by load_data.py; nothing here reads the CSV directly.
"""

import sqlite3
import warnings
from pathlib import Path

import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
ALPHA = 0.05


def connect(db_path=DB_PATH):
    if not Path(db_path).exists():
        raise FileNotFoundError(f"{db_path} not found - run `python load_data.py` first")
    return sqlite3.connect(db_path)


# --------------------------------------------------------------------------- #
# Part 2: relative frequency of each population in each sample
# --------------------------------------------------------------------------- #
FREQUENCY_SQL = """
WITH totals AS (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
)
SELECT
    cc.sample_id                          AS sample,
    t.total_count                         AS total_count,
    cp.name                               AS population,
    cc.count                              AS count,
    100.0 * cc.count / t.total_count      AS percentage
FROM cell_counts cc
JOIN totals t            ON t.sample_id = cc.sample_id
JOIN cell_populations cp ON cp.population_id = cc.population_id
ORDER BY cc.sample_id, cp.population_id
"""


def frequency_table(conn):
    """One row per (sample, population): sample, total_count, population, count, percentage."""
    return pd.read_sql_query(FREQUENCY_SQL, conn)


# --------------------------------------------------------------------------- #
# Part 3: responders vs non-responders (melanoma, miraclib, PBMC)
# --------------------------------------------------------------------------- #
RESPONSE_COHORT_SQL = f"""
WITH totals AS (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
)
SELECT
    s.sample_id                           AS sample,
    sub.subject_id                        AS subject,
    sub.project_id                        AS project,
    sub.response                          AS response,
    s.time_from_treatment_start           AS time_from_treatment_start,
    cp.name                               AS population,
    cc.count                              AS count,
    t.total_count                         AS total_count,
    100.0 * cc.count / t.total_count      AS percentage
FROM cell_counts cc
JOIN totals t            ON t.sample_id = cc.sample_id
JOIN cell_populations cp ON cp.population_id = cc.population_id
JOIN samples s           ON s.sample_id = cc.sample_id
JOIN subjects sub        ON sub.subject_id = s.subject_id
WHERE sub.condition = 'melanoma'
  AND sub.treatment = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND sub.response IN ('yes', 'no')
ORDER BY s.sample_id, cp.population_id
"""


def response_cohort(conn):
    """Part 2 frequencies restricted to melanoma / miraclib / PBMC samples with a known response."""
    df = pd.read_sql_query(RESPONSE_COHORT_SQL, conn)
    df["response_group"] = df["response"].map({"yes": "Responder", "no": "Non-responder"})
    return df


def _bh(pvalues):
    return multipletests(pvalues, method="fdr_bh")[1]


def compare_responders(cohort):
    """Per-population responder vs non-responder comparison.

    Each subject contributes three PBMC samples (days 0, 7 and 14), so samples
    are not independent. The primary test therefore is a linear mixed model

        percentage ~ response + C(time_from_treatment_start) + (1 | subject)

    i.e. a fixed effect for response, adjusting for timepoint, with a random
    intercept per subject. Benjamini-Hochberg FDR is applied across the five
    populations and a population is called significant if q < ALPHA.

    Two non-parametric Mann-Whitney U tests are reported as supporting evidence:
      * sample level (all samples, ignores the repeated measures);
      * subject level (one per-subject mean frequency per subject).
    """
    rows = []
    for pop, d in cohort.groupby("population", sort=False):
        r = d.loc[d.response == "yes", "percentage"]
        n = d.loc[d.response == "no", "percentage"]
        u, p = stats.mannwhitneyu(r, n, alternative="two-sided")

        subj = d.groupby(["subject", "response"], as_index=False)["percentage"].mean()
        rs = subj.loc[subj.response == "yes", "percentage"]
        ns = subj.loc[subj.response == "no", "percentage"]
        _, p_subj = stats.mannwhitneyu(rs, ns, alternative="two-sided")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = smf.mixedlm(
                "percentage ~ C(response, Treatment('no')) + C(time_from_treatment_start)",
                d, groups=d["subject"],
            ).fit(reml=True)
        term = "C(response, Treatment('no'))[T.yes]"

        rows.append({
            "population": pop,
            "n_responder_samples": len(r),
            "n_nonresponder_samples": len(n),
            "n_responder_subjects": len(rs),
            "n_nonresponder_subjects": len(ns),
            "responder_median_pct": r.median(),
            "nonresponder_median_pct": n.median(),
            "responder_mean_pct": r.mean(),
            "nonresponder_mean_pct": n.mean(),
            # mixed model: responder - non-responder difference, percentage points
            "effect_pct_points": fit.params[term],
            "ci95_low": fit.conf_int().loc[term, 0],
            "ci95_high": fit.conf_int().loc[term, 1],
            "p_value": fit.pvalues[term],
            # rank-biserial correlation: P(R > N) - P(R < N); range -1..1
            "rank_biserial": 2 * u / (len(r) * len(n)) - 1,
            "p_mannwhitney_samples": p,
            "p_mannwhitney_subjects": p_subj,
        })

    res = pd.DataFrame(rows)
    res["q_value"] = _bh(res["p_value"])
    res["q_mannwhitney_samples"] = _bh(res["p_mannwhitney_samples"])
    res["q_mannwhitney_subjects"] = _bh(res["p_mannwhitney_subjects"])
    res["significant"] = res["q_value"] < ALPHA
    return res


def compare_by_timepoint(cohort):
    """Mann-Whitney U at each timepoint separately (one sample per subject per timepoint)."""
    rows = []
    for (t, pop), d in cohort.groupby(["time_from_treatment_start", "population"], sort=False):
        r = d.loc[d.response == "yes", "percentage"]
        n = d.loc[d.response == "no", "percentage"]
        u, p = stats.mannwhitneyu(r, n, alternative="two-sided")
        rows.append({
            "time_from_treatment_start": t,
            "population": pop,
            "n_responder": len(r),
            "n_nonresponder": len(n),
            "responder_median_pct": r.median(),
            "nonresponder_median_pct": n.median(),
            "rank_biserial": 2 * u / (len(r) * len(n)) - 1,
            "p_value": p,
        })
    res = pd.DataFrame(rows).sort_values(["time_from_treatment_start", "population"])
    # FDR across all timepoint x population tests
    res["q_value"] = _bh(res["p_value"])
    res["significant"] = res["q_value"] < ALPHA
    return res.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Part 4: baseline melanoma / miraclib / PBMC subset
# --------------------------------------------------------------------------- #
BASELINE_FILTER = """
    FROM samples s
    JOIN subjects sub ON sub.subject_id = s.subject_id
    WHERE sub.condition = 'melanoma'
      AND sub.treatment = 'miraclib'
      AND s.sample_type = 'PBMC'
      AND s.time_from_treatment_start = 0
"""

BASELINE_SAMPLES_SQL = f"""
SELECT s.sample_id AS sample, sub.subject_id AS subject, sub.project_id AS project,
       sub.condition, sub.treatment, sub.response, sub.sex, sub.age,
       s.sample_type, s.time_from_treatment_start
{BASELINE_FILTER}
ORDER BY s.sample_id
"""

# LEFT JOIN from projects so projects with no matching samples show 0.
SAMPLES_PER_PROJECT_SQL = f"""
SELECT p.project_id AS project, COUNT(b.sample_id) AS n_samples
FROM projects p
LEFT JOIN (SELECT s.sample_id, sub.project_id {BASELINE_FILTER}) b
       ON b.project_id = p.project_id
GROUP BY p.project_id
ORDER BY p.project_id
"""

SUBJECTS_BY_RESPONSE_SQL = f"""
SELECT CASE sub.response WHEN 'yes' THEN 'responder'
                         WHEN 'no'  THEN 'non-responder'
                         ELSE 'unknown' END AS response,
       COUNT(DISTINCT sub.subject_id) AS n_subjects
{BASELINE_FILTER}
GROUP BY sub.response
ORDER BY sub.response DESC
"""

SUBJECTS_BY_SEX_SQL = f"""
SELECT CASE sub.sex WHEN 'M' THEN 'male' WHEN 'F' THEN 'female' ELSE 'unknown' END AS sex,
       COUNT(DISTINCT sub.subject_id) AS n_subjects
{BASELINE_FILTER}
GROUP BY sub.sex
ORDER BY sub.sex DESC
"""


def baseline_samples(conn):
    return pd.read_sql_query(BASELINE_SAMPLES_SQL, conn)


def baseline_samples_per_project(conn):
    return pd.read_sql_query(SAMPLES_PER_PROJECT_SQL, conn)


def baseline_subjects_by_response(conn):
    return pd.read_sql_query(SUBJECTS_BY_RESPONSE_SQL, conn)


def baseline_subjects_by_sex(conn):
    return pd.read_sql_query(SUBJECTS_BY_SEX_SQL, conn)


def format_p(p):
    return f"{p:.2e}" if p < 1e-3 else f"{p:.4f}"


if __name__ == "__main__":
    pd.set_option("display.width", 200, "display.max_columns", 30)
    with connect() as conn:
        cohort = response_cohort(conn)
        print(compare_responders(cohort).T)
        print(compare_by_timepoint(cohort))
        print(baseline_samples_per_project(conn))
        print(baseline_subjects_by_response(conn))
        print(baseline_subjects_by_sex(conn))
