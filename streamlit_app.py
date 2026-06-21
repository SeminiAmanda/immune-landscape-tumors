import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from lifelines import KaplanMeierFitter

# ── Page config ──
st.set_page_config(
    page_title="The Immune Landscape of Tumors",
    page_icon="🧬",
    layout="wide"
)

# ── Load data (cached so it only loads once) ──
@st.cache_data
def load_data():
    master = pd.read_csv("master_clinical_immune.csv", index_col=0)
    umap_df = pd.read_csv("umap_coords.csv", index_col=0)
    return master, umap_df

master, umap_df = load_data()

IMMUNE_COLS = ['CD8_T_cell', 'CD4_T_cell', 'Treg', 'NK_cell', 'B_cell',
               'M1_Macrophage', 'M2_Macrophage', 'Dendritic_cell',
               'Exhausted_T', 'PD_L1_axis']

CANCER_TYPES = sorted(master['cancer type abbreviation'].dropna().unique())

# ── Header ──
st.markdown(
    "<h1 style='text-align: center; color: #1a1a2e;'>The Immune Landscape of Tumors</h1>",
    unsafe_allow_html=True
)
st.markdown(
    "<p style='text-align: center; color: #666; font-size: 14px;'>"
    "A Multi-Omics Atlas of Cancer Microenvironment — CFDE Data Visualization Competition</p>",
    unsafe_allow_html=True
)

# ── Tabs ──
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🗺️ Immune Landscape",
    "🔥 Infiltration Heatmap",
    "🧬 Mutation Correlation",
    "⏳ Survival Explorer",
    "🔬 Gene Deep Dive"
])

# ════════════════════════════════════════════════
# TAB 1: Immune Landscape (UMAP)
# ════════════════════════════════════════════════
with tab1:
    st.info(
        "📖 **How to Read This View**  \n"
        "Each dot is one tumor sample. Tumors cluster by how much immune activity "
        "they show — 'hot' tumors (rich in immune cells) tend to respond better to "
        "immunotherapy than 'cold' tumors."
    )

    col1, col2 = st.columns(2)
    with col1:
        selected_cancers_1 = st.multiselect(
            "Filter by Cancer Type:", CANCER_TYPES, key="cancer_filter_1"
        )
    with col2:
        color_mode = st.radio(
            "Color By:",
            options=["Subtype_Immune_Model_Based", "cancer type abbreviation", "TMB"],
            format_func=lambda x: {
                "Subtype_Immune_Model_Based": "Immune Subtype",
                "cancer type abbreviation": "Cancer Type",
                "TMB": "Tumor Mutational Burden"
            }[x],
            horizontal=True
        )

    df1 = umap_df.copy()
    if selected_cancers_1:
        df1 = df1[df1['cancer type abbreviation'].isin(selected_cancers_1)]

    color_kwargs = {}
    if color_mode == "TMB":
        color_kwargs = dict(color_continuous_scale="Viridis",
                             range_color=[0, df1["TMB"].quantile(0.95)])

    fig1 = px.scatter(
        df1, x="UMAP1", y="UMAP2", color=color_mode,
        hover_data={"cancer type abbreviation": True, "TMB": ":.2f",
                    "UMAP1": False, "UMAP2": False},
        title="Immune Landscape of Tumors (UMAP Projection)",
        **color_kwargs
    )
    fig1.update_traces(marker=dict(size=5, opacity=0.7))
    fig1.update_layout(plot_bgcolor="white", paper_bgcolor="white", height=600,
                        legend_title_text=color_mode.replace("_", " "))
    st.plotly_chart(fig1, use_container_width=True)

# ════════════════════════════════════════════════
# TAB 2: Immune Infiltration Heatmap
# ════════════════════════════════════════════════
with tab2:
    st.info(
        "📖 **How to Read This View**  \n"
        "This heatmap shows average immune cell activity per cancer type. "
        "Brighter colors mean that immune cell type is more abundant in that cancer's "
        "tumors. Look for cancers with high CD8+ T cell / low Treg — these are good "
        "immunotherapy candidates."
    )

    avg_by_cancer = master.groupby("cancer type abbreviation")[IMMUNE_COLS].mean()
    z = (avg_by_cancer - avg_by_cancer.mean()) / avg_by_cancer.std()

    fig2 = go.Figure(data=go.Heatmap(
        z=z.values, x=IMMUNE_COLS, y=z.index,
        colorscale="RdBu_r", zmid=0,
        colorbar=dict(title="Z-score")
    ))
    fig2.update_layout(
        title="Immune Cell Infiltration by Cancer Type (Z-scored)",
        xaxis_title="Immune Cell Type", yaxis_title="Cancer Type",
        plot_bgcolor="white", paper_bgcolor="white", height=650
    )
    st.plotly_chart(fig2, use_container_width=True)

# ════════════════════════════════════════════════
# TAB 3: Mutation ↔ Immune Correlation
# ════════════════════════════════════════════════
with tab3:
    st.info(
        "📖 **How to Read This View**  \n"
        "Each dot is a tumor. The x-axis shows how many mutations it carries "
        "(Tumor Mutational Burden). The y-axis shows CD8+ T cell infiltration. "
        "Tumors with high TMB often have more immune cell infiltration — this is "
        "because more mutations create more 'foreign-looking' proteins the immune "
        "system can recognize."
    )

    selected_cancers_3 = st.multiselect(
        "Filter by Cancer Type:", CANCER_TYPES, key="cancer_filter_3"
    )

    df3 = master.dropna(subset=["TMB", "CD8_T_cell"]).copy()
    if selected_cancers_3:
        df3 = df3[df3["cancer type abbreviation"].isin(selected_cancers_3)]
    df3 = df3[df3["TMB"] < df3["TMB"].quantile(0.98)]

    fig3 = px.scatter(
        df3, x="TMB", y="CD8_T_cell", color="cancer type abbreviation",
        trendline="ols", trendline_scope="overall",
        hover_data={"Subtype_Immune_Model_Based": True},
        title="Tumor Mutational Burden vs CD8+ T Cell Infiltration",
        labels={"TMB": "Tumor Mutational Burden (mut/Mb)", "CD8_T_cell": "CD8+ T Cell Score"}
    )
    fig3.update_traces(marker=dict(size=5, opacity=0.6), selector=dict(mode="markers"))
    fig3.update_layout(plot_bgcolor="white", paper_bgcolor="white", height=600)
    st.plotly_chart(fig3, use_container_width=True)

# ════════════════════════════════════════════════
# TAB 4: Survival Explorer
# ════════════════════════════════════════════════
with tab4:
    st.info(
        "📖 **How to Read This View**  \n"
        "This curve shows the percentage of patients still alive over time. "
        "We split tumors into 'Hot' (high CD8+ T cell infiltration, top 50%) vs "
        "'Cold' (bottom 50%) and compare their survival. A curve that stays higher "
        "for longer means better survival."
    )

    selected_cancers_4 = st.multiselect(
        "Filter by Cancer Type:", CANCER_TYPES, key="cancer_filter_4"
    )

    df4 = master.dropna(subset=["OS", "OS.time", "CD8_T_cell"]).copy()
    if selected_cancers_4:
        df4 = df4[df4["cancer type abbreviation"].isin(selected_cancers_4)]

    median_cd8 = df4["CD8_T_cell"].median()
    df4["group"] = np.where(df4["CD8_T_cell"] >= median_cd8, "Hot (High CD8+)", "Cold (Low CD8+)")

    fig4 = go.Figure()
    for group, color in [("Hot (High CD8+)", "#d62728"), ("Cold (Low CD8+)", "#1f77b4")]:
        subset = df4[df4["group"] == group]
        kmf = KaplanMeierFitter()
        kmf.fit(subset["OS.time"], event_observed=subset["OS"], label=group)
        surv = kmf.survival_function_.reset_index()
        fig4.add_trace(go.Scatter(
            x=surv["timeline"], y=surv[group],
            mode="lines", name=f"{group} (n={len(subset)})",
            line=dict(color=color, width=2.5)
        ))

    fig4.update_layout(
        title="Survival: Immune-Hot vs Immune-Cold Tumors",
        xaxis_title="Days", yaxis_title="Survival Probability",
        plot_bgcolor="white", paper_bgcolor="white",
        yaxis=dict(range=[0, 1.05]), height=600
    )
    st.plotly_chart(fig4, use_container_width=True)

# ════════════════════════════════════════════════
# TAB 5: Gene Expression Deep Dive
# ════════════════════════════════════════════════
with tab5:
    st.info(
        "📖 **How to Read This View**  \n"
        "Each box shows the spread of expression for one immune cell signature "
        "across immune subtypes. Wider separation between boxes means that signature "
        "strongly distinguishes hot vs cold tumors."
    )

    selected_signature = st.selectbox(
        "Select Immune Signature:",
        options=IMMUNE_COLS,
        format_func=lambda c: c.replace("_", " "),
        index=0
    )

    df5 = master.dropna(subset=[selected_signature, "Subtype_Immune_Model_Based"]).copy()

    fig5 = px.violin(
        df5, x="Subtype_Immune_Model_Based", y=selected_signature,
        color="Subtype_Immune_Model_Based", box=True, points=False,
        title=f'{selected_signature.replace("_", " ")} Score by Immune Subtype'
    )
    fig5.update_layout(
        plot_bgcolor="white", paper_bgcolor="white",
        showlegend=False, xaxis_title="Immune Subtype",
        yaxis_title=f'{selected_signature.replace("_", " ")} Score',
        xaxis_tickangle=-20, height=600
    )
    st.plotly_chart(fig5, use_container_width=True)
