"""Streamlit dashboard for the miraclib analysis.

Start it with `make dashboard`. Needs cell_counts.db from `make pipeline`.
"""

import plotly.graph_objects as go
import streamlit as st

import analysis

st.set_page_config(page_title="Miraclib Immune Profiling", layout="wide")


# streamlit runs sessions on different threads, so turn off sqlite's thread check
# (fine since we only read)
@st.cache_resource
def get_conn():
    return analysis.connect(check_same_thread=False)


@st.cache_data
def load_frequency():
    return analysis.frequency_table(get_conn())


@st.cache_data
def load_responders():
    return analysis.responder_frequencies(get_conn())


@st.cache_data
def load_subset():
    return analysis.baseline_subset(get_conn())


try:
    get_conn()
except FileNotFoundError as e:
    st.error(f"{e} Run `make pipeline` first.")
    st.stop()

st.title("Miraclib trial: immune cell populations")
tab2, tab3, tab4 = st.tabs(["Part 2: Data overview", "Part 3: Responders", "Part 4: Baseline subset"])

with tab2:
    st.subheader("Relative frequency of each population per sample")
    freq = load_frequency()

    c1, c2, c3 = st.columns(3)
    c1.metric("Samples", f"{freq['sample'].nunique():,}")
    c2.metric("Populations", freq["population"].nunique())
    c3.metric("Rows", f"{len(freq):,}")

    f1, f2 = st.columns([2, 1])
    search = f1.text_input("Search sample id", "")
    pops = f2.multiselect("Populations", analysis.POPULATIONS, default=analysis.POPULATIONS)

    view = freq[freq["population"].isin(pops)]
    if search:
        view = view[view["sample"].str.contains(search.strip(), case=False)]
    st.dataframe(view, hide_index=True, width="stretch",
                 column_config={"percentage": st.column_config.NumberColumn(format="%.2f")})
    st.download_button("Download CSV", freq.to_csv(index=False), "frequency_table.csv", "text/csv")

with tab3:
    st.subheader("Melanoma, miraclib, PBMC: responders vs non-responders")
    resp = load_responders()
    results = analysis.responder_statistics(resp)

    labels = {
        "all_samples": "All samples",
        "subject_mean": "Mean per subject",
        "baseline": "Baseline only (day 0)",
    }
    mode = st.radio("Analysis", list(results), format_func=labels.get, horizontal=True)

    if mode == "subject_mean":
        plot_df = resp.groupby(["subject_id", "population", "group"], as_index=False)["percentage"].mean()
    elif mode == "baseline":
        plot_df = resp[resp["time_from_treatment_start"] == 0]
    else:
        plot_df = resp

    fig = go.Figure()
    for group, color in analysis.COLORS.items():
        d = plot_df[plot_df["group"] == group]
        fig.add_trace(go.Box(x=d["population"], y=d["percentage"], name=group, marker_color=color))
    fig.update_layout(boxmode="group", height=480, margin=dict(t=30, b=40),
                      xaxis_title="Cell population", yaxis_title="Relative frequency (%)")
    st.plotly_chart(fig, width="stretch")

    st.markdown(
        "Two-sided Mann-Whitney U test for each population, with p-values adjusted for the "
        f"5 comparisons using Benjamini-Hochberg (q). Threshold is {analysis.ALPHA}. "
        "rank_biserial is the effect size (positive = higher in responders)."
    )
    table = results[mode]
    num_cols = ["median_responders", "median_non_responders", "median_diff",
                "rank_biserial", "p_value", "q_value"]
    st.dataframe(table, hide_index=True, width="stretch",
                 column_config={c: st.column_config.NumberColumn(format="%.4f") for c in num_cols})

    sig = table.loc[table["significant"], "population"].tolist()
    nominal = table.loc[table["significant_nominal"], "population"].tolist()
    if sig:
        st.success(f"Significant after correction (q < {analysis.ALPHA}): {', '.join(sig)}")
    else:
        st.warning(f"Nothing is significant after correction (q < {analysis.ALPHA}).")
    if nominal and set(nominal) != set(sig):
        st.info(f"p < {analysis.ALPHA} before correction: {', '.join(nominal)}")

with tab4:
    st.subheader("Melanoma, miraclib, PBMC samples at baseline (day 0)")
    subset = load_subset()
    st.metric("Baseline samples", len(subset["samples"]))

    c1, c2, c3 = st.columns(3)
    c1.markdown("**Samples per project**")
    c1.dataframe(subset["samples_per_project"], hide_index=True, width="stretch")
    c2.markdown("**Subjects by response**")
    c2.dataframe(subset["subjects_by_response"], hide_index=True, width="stretch")
    c3.markdown("**Subjects by sex**")
    c3.dataframe(subset["subjects_by_sex"], hide_index=True, width="stretch")

    with st.expander("Show all baseline samples"):
        st.dataframe(subset["samples"], hide_index=True, width="stretch")

    st.divider()
    avg_b = analysis.avg_b_cells_melanoma_male_responders_t0(get_conn())
    st.metric("Avg B cells: melanoma males, responders, day 0 (all sample types/treatments)", f"{avg_b:.2f}")
