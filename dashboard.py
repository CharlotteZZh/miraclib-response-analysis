"""Interactive dashboard for the Part 2-4 results.

Run with:  make dashboard   (or: streamlit run dashboard.py)
Reads everything from cell_counts.db created by `make pipeline`.
"""

import plotly.express as px
import streamlit as st

import analysis

GROUP_ORDER = ["Non-responder", "Responder"]
GROUP_COLORS = {"Non-responder": "#C2410C", "Responder": "#2563EB"}

st.set_page_config(page_title="Miraclib Response Analysis", layout="wide")


@st.cache_data
def load():
    with analysis.connect() as conn:
        cohort = analysis.response_cohort(conn)
        return {
            "freq": analysis.frequency_table(conn),
            "cohort": cohort,
            "stats": analysis.compare_responders(cohort),
            "by_time": analysis.compare_by_timepoint(cohort),
            "baseline": analysis.baseline_samples(conn),
            "per_project": analysis.baseline_samples_per_project(conn),
            "by_response": analysis.baseline_subjects_by_response(conn),
            "by_sex": analysis.baseline_subjects_by_sex(conn),
        }


try:
    data = load()
except FileNotFoundError as e:
    st.error(f"{e}. Run `make pipeline` first.")
    st.stop()

st.title("Miraclib: immune cell populations and treatment response")
st.caption("Loblaw Bio clinical trial - cell-count.csv loaded into SQLite (cell_counts.db)")

tab2, tab3, tab4 = st.tabs([
    "Part 2 - Frequency overview",
    "Part 3 - Responders vs non-responders",
    "Part 4 - Baseline subset",
])

# --------------------------------------------------------------------------- #
with tab2:
    freq = data["freq"]
    st.subheader("Relative frequency of each cell population per sample")
    c1, c2, c3 = st.columns(3)
    c1.metric("Samples", f"{freq['sample'].nunique():,}")
    c2.metric("Populations", freq["population"].nunique())
    c3.metric("Rows", f"{len(freq):,}")

    f1, f2 = st.columns([2, 1])
    sample_query = f1.text_input("Filter by sample id (substring)", "")
    pops = f2.multiselect("Populations", freq["population"].unique().tolist(),
                          default=freq["population"].unique().tolist())
    view = freq[freq["population"].isin(pops)]
    if sample_query:
        view = view[view["sample"].str.contains(sample_query.strip(), case=False)]
    st.dataframe(
        view, hide_index=True, width="stretch", height=420,
        column_config={"percentage": st.column_config.NumberColumn(format="%.2f %%")},
    )
    st.download_button("Download full table (CSV)", freq.to_csv(index=False),
                       "part2_frequency_table.csv", "text/csv")

    st.markdown("**Distribution of relative frequencies across all samples**")
    fig = px.box(freq, x="population", y="percentage", points=False,
                 labels={"percentage": "Relative frequency (%)", "population": ""})
    st.plotly_chart(fig, width="stretch")

# --------------------------------------------------------------------------- #
with tab3:
    cohort, stats_df = data["cohort"], data["stats"]
    st.subheader("Melanoma patients on miraclib, PBMC samples only")
    n_r = cohort.loc[cohort.response == "yes", "subject"].nunique()
    n_n = cohort.loc[cohort.response == "no", "subject"].nunique()
    c1, c2, c3 = st.columns(3)
    c1.metric("Samples", f"{cohort['sample'].nunique():,}")
    c2.metric("Responder subjects", n_r)
    c3.metric("Non-responder subjects", n_n)

    timepoints = sorted(cohort["time_from_treatment_start"].unique())
    tsel = st.multiselect("Timepoints (days from treatment start) shown in the plot",
                          timepoints, default=timepoints)
    plot_df = cohort[cohort["time_from_treatment_start"].isin(tsel)]
    fig = px.box(
        plot_df, x="population", y="percentage", color="response_group",
        category_orders={"response_group": GROUP_ORDER,
                         "population": stats_df["population"].tolist()},
        color_discrete_map=GROUP_COLORS, points=False,
        labels={"percentage": "Relative frequency (%)", "population": "",
                "response_group": ""},
    )
    fig.update_layout(boxmode="group", legend={"orientation": "h", "y": 1.08})
    st.plotly_chart(fig, width="stretch")

    sig = stats_df.loc[stats_df.significant, "population"].tolist()
    if sig:
        st.success(
            "Significant difference (mixed model, BH-adjusted q < 0.05): **"
            + ", ".join(sig) + "**"
        )
    else:
        st.info("No population differs significantly after FDR correction.")

    st.markdown("**Statistical results (all timepoints)**")
    show = stats_df[[
        "population", "responder_mean_pct", "nonresponder_mean_pct",
        "effect_pct_points", "ci95_low", "ci95_high", "p_value", "q_value",
        "p_mannwhitney_samples", "p_mannwhitney_subjects", "rank_biserial", "significant",
    ]]
    st.dataframe(show.style.format(precision=4), hide_index=True, width="stretch")
    st.markdown(
        """
- **Primary test - linear mixed model:** `percentage ~ response + timepoint + (1 | subject)`.
  Each subject has three PBMC samples (day 0, 7, 14), so a random intercept per
  subject accounts for the repeated measures. `effect_pct_points` is the responder
  minus non-responder difference in percentage points, with its 95% CI.
  `q_value` is the Benjamini-Hochberg-adjusted p-value across the 5 populations.
- **Supporting tests - Mann-Whitney U (two-sided):** on all samples, and on
  per-subject mean frequencies (one value per subject). `rank_biserial` is the
  effect size (-1 to 1).
"""
    )

    st.markdown("**Per-timepoint Mann-Whitney U tests** (BH-adjusted across all 15 tests)")
    st.dataframe(data["by_time"].style.format(precision=4), hide_index=True,
                 width="stretch")

# --------------------------------------------------------------------------- #
with tab4:
    st.subheader("Melanoma PBMC samples at baseline (day 0) from miraclib-treated patients")
    base = data["baseline"]
    c1, c2, c3 = st.columns(3)
    c1.metric("Samples", len(base))
    c2.metric("Subjects", base["subject"].nunique())
    c3.metric("Projects with samples", int((data["per_project"]["n_samples"] > 0).sum()))

    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**Samples per project**")
        st.dataframe(data["per_project"], hide_index=True, width="stretch")
        st.plotly_chart(px.bar(data["per_project"], x="project", y="n_samples", height=260),
                        width="stretch")
    with c2:
        st.markdown("**Subjects by response**")
        st.dataframe(data["by_response"], hide_index=True, width="stretch")
        st.plotly_chart(px.bar(data["by_response"], x="response", y="n_subjects", height=260),
                        width="stretch")
    with c3:
        st.markdown("**Subjects by sex**")
        st.dataframe(data["by_sex"], hide_index=True, width="stretch")
        st.plotly_chart(px.bar(data["by_sex"], x="sex", y="n_subjects", height=260),
                        width="stretch")

    with st.expander("All baseline samples"):
        st.dataframe(base, hide_index=True, width="stretch")
