#!/usr/bin/env python3
"""Run the full pipeline end to end: build the database (Part 1), then compute
and save the Part 2-4 results to outputs/.

Usage:
    python run_pipeline.py
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import analysis  # noqa: E402
import load_data  # noqa: E402

OUT = analysis.ROOT / "outputs"
GROUP_ORDER = ["Non-responder", "Responder"]
GROUP_COLORS = {"Non-responder": "#C2410C", "Responder": "#2563EB"}


def boxplot(cohort, stats_df, path):
    pops = list(stats_df["population"])
    fig, axes = plt.subplots(1, len(pops), figsize=(3.2 * len(pops), 4.6))
    for ax, pop in zip(axes, pops):
        d = cohort[cohort.population == pop]
        data = [d.loc[d.response_group == g, "percentage"] for g in GROUP_ORDER]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=True,
                        flierprops={"markersize": 2, "alpha": 0.4})
        for patch, g in zip(bp["boxes"], GROUP_ORDER):
            patch.set_facecolor(GROUP_COLORS[g])
            patch.set_alpha(0.55)
        for median in bp["medians"]:
            median.set_color("black")
        row = stats_df.set_index("population").loc[pop]
        star = " *" if row["significant"] else ""
        ax.set_title(f"{pop}{star}\nq = {analysis.format_p(row['q_value'])}", fontsize=10)
        ax.set_xticks([1, 2], GROUP_ORDER, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Relative frequency (% of total cells)")
    fig.suptitle(
        "Melanoma, miraclib, PBMC: responders vs non-responders "
        "(q = BH-adjusted mixed-model p; * q < 0.05)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    # Part 1
    load_data.load()
    OUT.mkdir(exist_ok=True)

    with analysis.connect() as conn:
        # Part 2
        freq = analysis.frequency_table(conn)
        freq.to_csv(OUT / "part2_frequency_table.csv", index=False)
        print(f"Part 2: {len(freq)} rows -> outputs/part2_frequency_table.csv")

        # Part 3
        cohort = analysis.response_cohort(conn)
        stats_df = analysis.compare_responders(cohort)
        stats_df.to_csv(OUT / "part3_responder_stats.csv", index=False)
        analysis.compare_by_timepoint(cohort).to_csv(
            OUT / "part3_responder_stats_by_timepoint.csv", index=False)
        boxplot(cohort, stats_df, OUT / "part3_boxplot.png")
        sig = stats_df.loc[stats_df.significant, "population"].tolist()
        print("Part 3: significant populations (mixed model, BH q < 0.05):",
              ", ".join(sig) if sig else "none")
        print(stats_df[["population", "effect_pct_points", "p_value", "q_value",
                        "p_mannwhitney_samples", "p_mannwhitney_subjects"]]
              .to_string(index=False, float_format=lambda x: f"{x:.4g}"))

        # Part 4
        analysis.baseline_samples(conn).to_csv(OUT / "part4_baseline_samples.csv", index=False)
        for name, fn in [
            ("part4_samples_per_project", analysis.baseline_samples_per_project),
            ("part4_subjects_by_response", analysis.baseline_subjects_by_response),
            ("part4_subjects_by_sex", analysis.baseline_subjects_by_sex),
        ]:
            df = fn(conn)
            df.to_csv(OUT / f"{name}.csv", index=False)
            print(f"{name}:\n{df.to_string(index=False)}")

    print("Pipeline complete; results in outputs/")


if __name__ == "__main__":
    main()
