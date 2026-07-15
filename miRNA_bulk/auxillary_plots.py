import warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve
from constants import *
from matplotlib.colors import LinearSegmentedColormap


FINAL_CHOICES = [
        'Tissue specific miR',
        'Cancer-tissue-specific miR',
        'Normal-tissue-specific miR',
        'Non-specific miR'
    ]

FINAL_COLORS = {
    'Tissue specific miR': '#2ca02c',
    'Cancer-tissue-specific miR': '#9467bd',
    'Normal-tissue-specific miR': '#B8860B',
    'Non-specific miR': '#7f7f7f'
}

DYSREGULATED_CLASS_COLORS = {
    "NOT_SIGNIFICANT": "gray",
    "UPREGULATED": "firebrick",
    "DOWNREGULATED": "navy",
}

STN_CLASS_COLORS = {
    "NOT_SIGNIFICANT": "gray",
    "UPREGULATED": "firebrick",
    "DOWNREGULATED": "navy",
}

BINOM_COLORS = {
    "UP": "firebrick",
    "DOWN": "navy",
    "NS": "gray",
}

DIRECTION_COLORS = {
    "TCGA_UP": "firebrick",
    "GTEx_UP": "navy",
    "NS": "gray",
    "NA": "#cccccc",
}

neutral_cutoff = 0.25
vmax = 20
w = neutral_cutoff / vmax
_DELTA_GAIN_COLOR = "#9467BD"  
_DELTA_LOSS_COLOR = "#DC143C" 

colors = [
    (0.00, "#08306B"),
    (0.22, "#4292C6"),

    (0.5 - w, "white"), #f2e8d5
    (0.5 + w, "white"), #f2e8d5

    (0.80, "#FC9272"),
    (1.00, "#99000D")
]

custom_cmap = LinearSegmentedColormap.from_list(
    "de_map",
    colors
)

def _set_specificity_scatter_axes(ax, *, pad: float = SPECIFICITY_AXIS_PAD, equal_aspect: bool = True) -> None:
    """Pad τ axes so markers at x/y=0 are not clipped at the plot border."""
    lo, hi = -pad, 1.0 + pad
    ax.set_autoscale_on(False)
    ax.margins(x=0, y=0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    for coll in ax.collections:
        coll.set_clip_on(False)


def plot_correlation(df_in, save_path=None):
    delta_df = df_in.copy()
    x = delta_df["log2FC_tcga_gtex"]
    y = delta_df[f"Delta_{SPECIFICITY_METRIC}"]

    p, _ = pearsonr(x, y)
    p_sp, _ = spearmanr(x, y)

    plt.figure(figsize=(6,5))

    sns.scatterplot(
        data=delta_df,
        x="log2FC_tcga_gtex",
        y=f"Delta_{SPECIFICITY_METRIC}",
        alpha=0.5,
        edgecolor="black",
        linewidth=0.3
    )

    # LOWESS trend
    sns.regplot(
        data=delta_df,
        x="log2FC_tcga_gtex",
        y=f"Delta_{SPECIFICITY_METRIC}",
        scatter=False,
        lowess=True
    )

    # zero lines
    plt.axhline(0, linestyle="--", linewidth=1)
    plt.axvline(0, linestyle="--", linewidth=1)

    plt.xlabel("log2FC (TCGA - GTEx)")
    plt.ylabel(f"Δ{SPECIFICITY_METRIC} (TCGA - GTEx)")

    plt.text(
        0.65, 0.95,
        f"Pearson r = {p:.2f}\nSpearman r = {p_sp:.2f}",
        transform=plt.gca().transAxes,
        verticalalignment="top"
    )

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()

def plot_specificity_distribution(df: pd.DataFrame, x: str, title: str=f'{SPECIFICITY_METRIC} Distribution', save_path: str = None):
    
    sns.set_theme(style="whitegrid", palette="muted")
    plt.figure(figsize=(9, 5))

    ax = sns.kdeplot(
        data=df,
        x=x,
        fill=True,            
        color="#2c7fb8",      
        alpha=0.6,            
        linewidth=2,
        bw_adjust=0.5         
    )

    sns.rugplot(
        data=df,
        x=x,
        color="#2c7fb8",
        alpha=0.1,           
        height=0.05
    )

    plt.title(title, fontsize=14, pad=15, fontweight='bold')
    plt.xlabel(f"{SPECIFICITY_METRIC} score", fontsize=12, labelpad=10)
    plt.ylabel("Density", fontsize=12, labelpad=10)

    plt.xlim(-0.05, 1.05)
    sns.despine(left=True, bottom=True)
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", facecolor=plt.gcf().get_facecolor())
    plt.tight_layout()
    plt.show()



def counts_for_pie(series: pd.Series, top_n: int = TOP_N):
    s = series.fillna(MISSING).astype(str)
    vc = s.value_counts()

    if len(vc) <= top_n:
        return vc.index.tolist(), vc.values.astype(float)

    top = vc.iloc[:top_n]
    rest = vc.iloc[top_n:].sum()

    labels = top.index.tolist() + [OTHER]
    sizes = np.append(top.values.astype(float), float(rest))

    return labels, sizes


def _cmap(name):
    try:
        return plt.colormaps[name]
    except:
        return plt.cm.get_cmap(name)

def plot_pie_single(
    data: pd.DataFrame,
    col: str,
    title_prefix: str = "Samples",
    save_path: str = None,
):
    
    fig, ax = plt.subplots(figsize=(10, 8), dpi=120)
    fig.patch.set_facecolor("#fafafa")

    cmap = _cmap("tab20b")

    if col not in data.columns:
        raise ValueError(f"{col} not found in dataframe")

    labels, sizes = counts_for_pie(data[col], TOP_N)

    n = len(labels)
    colors = [cmap(i / max(n-1, 1)) for i in range(n)]

    for i, lab in enumerate(labels):
        if lab in (OTHER, MISSING):
            colors[i] = (0.82, 0.82, 0.86, 1.0)

    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=None,
        autopct=lambda p: f"{p:.1f}%" if p >= 2 else "",
        pctdistance=0.72,
        colors=colors,
        wedgeprops=dict(
            width=0.45,
            edgecolor="black",
            linewidth=0.8
        ),
        startangle=90,
    )

    for t in autotexts:
        t.set_fontsize(9)
        t.set_color("#222")

    ax.set_title(
        f"{title_prefix}: {col}",
        fontsize=13,
        fontweight="600",
        pad=16
    )

    leg_labels = [
        f"{lab} (n={int(sz)})"
        for lab, sz in zip(labels, sizes)
    ]

    ax.legend(
        wedges,
        leg_labels,
        title="Category",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=10,
        title_fontsize=10,
        frameon=False,
    )

    plt.suptitle(
        f"Top {TOP_N} + {OTHER} (by sample count)",
        fontsize=14
    )

    plt.tight_layout()

    if save_path:
        plt.savefig(
            save_path,
            bbox_inches="tight",
            facecolor=fig.get_facecolor()
        )

    plt.show()


def _get_val(r, col, is_padj=False):
    v = r.get(col)
    if pd.isna(v) or v == "–":
        return "NA"
    if isinstance(v, (int, float, np.floating)):
        return f"{v:.4e}" if is_padj else f"{v:.4f}"
    try:
        v = float(v)
        return f"{v:.4e}" if is_padj else f"{v:.4f}"
    except (TypeError, ValueError):
        return str(v)


def _fmt_cpm(v):
    if pd.isna(v) or v == "–":
        return "NA"
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if abs(v) >= 100:
        return f"{v:.1f}"
    if abs(v) >= 0.01:
        return f"{v:.4f}"
    return f"{v:.4e}"


def _mirna_hover_text(r, *, include_delta: bool = True) -> str:
    lines = [f"<b>{r['mature_name']}</b>", ""]
    lines.append("<b>Specificity:</b>")
    ann_tcga = f" ({r['annotation_tcga']})" if pd.notna(r.get("annotation_tcga")) else ""
    lines.append(f"{SPECIFICITY_METRIC} TCGA: {r[f'{SPECIFICITY_METRIC}_tcga']:.4f}{ann_tcga}")
    ann_gtex = f" ({r['annotation_gtex']})" if pd.notna(r.get("annotation_gtex")) else ""
    lines.append(f"{SPECIFICITY_METRIC} GTEx: {r[f'{SPECIFICITY_METRIC}_gtex']:.4f}{ann_gtex}")
    if include_delta and "spec_delta" in r.index and pd.notna(r.get("spec_delta")):
        lines.append(f"Δ {SPECIFICITY_METRIC} (TCGA−GTEx): {r['spec_delta']:+.4f}")
    lines.append(f"Tissue specificity (TCGA): {r.get('tcga_tissues', 'NA')}")
    lines.append(f"Tissue specificity (GTEx): {r.get('gtex_tissues', 'NA')}")
    lines.append("-----------------------")
    lines.append("<b>GTEx vs TCGA (pan-tissue CPM):</b>")
    lines.append(f"mean CPM GTEx: {_fmt_cpm(r.get('mean_cpm_gtex'))}")
    lines.append(f"max CPM GTEx: {_fmt_cpm(r.get('max_cpm_gtex'))}")
    lines.append(f"mean CPM TCGA: {_fmt_cpm(r.get('mean_cpm_tcga'))}")
    lines.append(f"max CPM TCGA: {_fmt_cpm(r.get('max_cpm_tcga'))}")
    lines.append(f"log2FC (TCGA/GTEx): {_get_val(r, 'log2FC_tcga_gtex')}")
    lines.append(f"direction: {r.get('direction', 'NA')}")
    lines.append(f"padj (MW): {_get_val(r, 'padj_tcga_gtex', True)}")
    lines.append("-----------------------")
    lines.append("<b>DE cancer vs normal:</b>")
    lines.append(f"n_projects_up: {r.get('n_projects_up', 'NA')}")
    lines.append(f"n_projects_down: {r.get('n_projects_down', 'NA')}")
    lines.append(f"median_log2FC: {_get_val(r, 'binom_median_log2fc')}")
    lines.append(f"median padj: {_get_val(r, 'binom_median_padj', True)}")
    lines.append(f"Binomial result: {r.get('Binomial_result', 'NA')}")
    lines.append("-----------------------")
    dysreg_class = r.get("Dysregulated_class", "NA")
    dysreg_class = ' '.join(dysreg_class.split('_'))
    lines.append(f"Dysregulation class: <b>{dysreg_class}</b>")
    return "<br>".join(lines)




def _ordered_categories(seen: pd.Series, preferred: list[str]) -> list[str]:
    cats = [c for c in preferred if (seen == c).any()]
    for c in sorted(seen.unique()):
        if c not in cats:
            cats.append(c)
    return cats


def _scatter_traces_by_column(
    df: pd.DataFrame,
    col: str,
    categories: list[str],
    palette: dict,
    hover: list[str],
    *,
    visible: bool,
    mk: dict,
) -> list:
    seen = df[col].astype(str).fillna("NS")
    traces = []
    for cat in categories:
        sub = df[seen == cat]
        if sub.empty:
            continue
        traces.append(go.Scatter(
            x=sub[f"{SPECIFICITY_METRIC}_gtex"],
            y=sub[f"{SPECIFICITY_METRIC}_tcga"],
            mode="markers",
            name=cat,
            visible=visible,
            marker={**mk, "color": palette.get(cat, "#888888")},
            text=[hover[i] for i in sub.index],
            hoverinfo="text",
        ))
    return traces


def _log2fc_color_limit(vals, color_pct: float = 95) -> float:
    """Symmetric ±limit for log2FC color scale (matches interactive scatter)."""
    s = pd.to_numeric(vals, errors="coerce").dropna()
    if s.empty:
        return 2.0
    lim = float(np.percentile(np.abs(s), color_pct))
    return max(lim, 0.25)


def _scatter_trace_log2fc(
    df: pd.DataFrame,
    col: str,
    hover: list[str],
    *,
    visible: bool,
    mk: dict,
    color_pct: float = 95,
) -> go.Scatter:
    """log2FC gradient; color scale clipped at ±percentile(|log2FC|) to reduce outlier compression."""
    if col not in df.columns:
        raise ValueError(f"Column {col!r} not found in merged dataframe.")

    lim = _log2fc_color_limit(df[col], color_pct)
    cmin, cmax = -lim, lim

    return go.Scatter(
        x=df[f"{SPECIFICITY_METRIC}_gtex"],
        y=df[f"{SPECIFICITY_METRIC}_tcga"],
        mode="markers",
        name="log2FC (TCGA/GTEx)",
        visible=visible,
        marker={
            **mk,
            "color": df[col],
            "colorscale": "Viridis",
            "cmin": cmin,
            "cmax": cmax,
            "cmid": 0,
            "colorbar": dict(
                title=f"log2FC<br>(TCGA/GTEx)<br>[±{lim:.1f} @ p{color_pct:.0f}]",
                x=1.02,
                len=0.75,
                thickness=18,
            ),
            "showscale": True,
        },
        text=hover,
        hoverinfo="text",
        showlegend=False,
    )


def plot_specificity_scatter_interactive(merged, save_path=None):
    plot_df = merged.dropna(
        subset=[f"{SPECIFICITY_METRIC}_tcga", f"{SPECIFICITY_METRIC}_gtex", "Specificity_class", f'Delta_{SPECIFICITY_METRIC}']
    ).copy().reset_index(drop=True)
    plot_df.rename(columns={f'Delta_{SPECIFICITY_METRIC}': 'spec_delta'}, inplace=True)

    hover = [_mirna_hover_text(plot_df.loc[i]) for i in range(len(plot_df))]
    mk = dict(size=14, opacity=0.85, line=dict(width=0.3, color="black"))

    views = [
        dict(
            label="Specificity class",
            col="Specificity_class",
            preferred=FINAL_CHOICES,
            palette=FINAL_COLORS,
            title=f"TCGA vs GTEx {SPECIFICITY_METRIC} — specificity class (n={len(plot_df)})",
        ),
        dict(
            label="Binomial (TCGA DE)",
            col="Binomial_result",
            preferred=["UP", "DOWN", "NS"],
            palette=BINOM_COLORS,
            title=f"TCGA vs GTEx {SPECIFICITY_METRIC} — Binomial_result (n={len(plot_df)})",
        ),
        dict(
            label="GTEx vs TCGA shift",
            col="direction",
            preferred=["TCGA_UP", "GTEx_UP", "NS"],
            palette=DIRECTION_COLORS,
            title=f"TCGA vs GTEx {SPECIFICITY_METRIC} — direction (n={len(plot_df)})",
        ),
    ]

    dysregulated_view = dict(
        label="Dysregulated",
        col="Dysregulated_class",
        preferred=["UPREGULATED", "DOWNREGULATED", "NOT_SIGNIFICANT"],
        palette=STN_CLASS_COLORS,
        title=(
            f"TCGA vs GTEx {SPECIFICITY_METRIC} — Dysregulated_class "
            f"(Binomial ∩ direction, n={len(plot_df)})"
        ),
    )

    log2fc_col = "log2FC_tcga_gtex"
    plot_df_lfc = plot_df[~plot_df["direction"].isin(["NS", "NA"])].reset_index(drop=True)
    hover_lfc = [_mirna_hover_text(plot_df_lfc.loc[i]) for i in range(len(plot_df_lfc))]

    log2fc_title = (
        f"TCGA vs GTEx {SPECIFICITY_METRIC} — log2FC TCGA/GTEx "
        f"(direction ≠ NS, n={len(plot_df_lfc)})"
    )
    view_titles = [v["title"] for v in views] + [log2fc_title, dysregulated_view["title"]]
    view_labels = [v["label"] for v in views] + ["log2FC (TCGA/GTEx)", dysregulated_view["label"]]

    trace_groups = []
    for i, view in enumerate(views):
        seen = plot_df[view["col"]].astype(str).fillna("NS")
        cats = _ordered_categories(seen, view["preferred"])
        trace_groups.append(
            _scatter_traces_by_column(
                plot_df, view["col"], cats, view["palette"], hover,
                visible=(i == 0), mk=mk,
            )
        )

    trace_groups.append([
        _scatter_trace_log2fc(
            plot_df_lfc, log2fc_col, hover_lfc, visible=False, mk=mk,
        )
    ])

    seen_dys = plot_df[dysregulated_view["col"]].astype(str).fillna("NOT_SIGNIFICANT")
    dys_cats = _ordered_categories(seen_dys, dysregulated_view["preferred"])
    trace_groups.append(
        _scatter_traces_by_column(
            plot_df, dysregulated_view["col"], dys_cats, dysregulated_view["palette"],
            hover, visible=False, mk=mk,
        )
    )

    traces = [t for group in trace_groups for t in group]
    group_sizes = [len(g) for g in trace_groups]

    fig = go.Figure(data=traces)

    _specificity_lo = -SPECIFICITY_AXIS_PAD
    _specificity_hi = 1.0 + SPECIFICITY_AXIS_PAD

    fig.add_shape(
        type="line", x0=_specificity_lo, y0=_specificity_lo, x1=_specificity_hi, y1=_specificity_hi,
        line=dict(color="rgba(0,0,0,0.4)", dash="dash"),
    )

    fig.update_xaxes(
        range=[_specificity_lo, _specificity_hi], title=f"{SPECIFICITY_METRIC} GTEx", autorange=False,
        constrain="domain",
        tickmode="array", tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
    )
    fig.update_yaxes(
        range=[_specificity_lo, _specificity_hi], title=f"{SPECIFICITY_METRIC} TCGA",
        scaleanchor="x", scaleratio=1, autorange=False, constrain="domain",
        tickmode="array", tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
    )

    def _visibility_mask(active_idx: int) -> list[bool]:
        vis = []
        for i, n in enumerate(group_sizes):
            vis.extend([i == active_idx] * n)
        return vis

    fig.update_layout(
        title=view_titles[0],
        template="plotly_white",
        height=850,
        autosize=True,
        margin=dict(l=70, r=280, t=120, b=70),
        legend=dict(
            orientation="v", x=1.02, y=1,
            xanchor="left", yanchor="top",
            bordercolor="LightGrey", borderwidth=1,
        ),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0, y=1.08,
                xanchor="left", yanchor="bottom",
                buttons=[
                    dict(
                        label=label,
                        method="update",
                        args=[
                            {"visible": _visibility_mask(i)},
                            {"title.text": view_titles[i]},
                        ],
                    )
                    for i, label in enumerate(view_labels)
                ],
            )
        ],
    )

    if save_path:
        fig.write_html(save_path, include_plotlyjs="cdn")

    fig.show()






def plot_stn_donut(
    df,
    class_col="Dysregulated_class",
    save_pdf=None
):

    counts = (
        df[class_col]
        .value_counts()
        .reindex([
            "UPREGULATED",
            "DOWNREGULATED",
            "NOT_SIGNIFICANT"
        ])
        .fillna(0)
    )

    labels = [
        "Upregulated",
        "Downregulated",
        "Non-significant"
    ]

    colors = [
        STN_CLASS_COLORS[k]
        for k in counts.index
    ]

    total = counts.sum()

    fig,ax = plt.subplots(
        figsize=(7,7),
        dpi=300
    )

    wedges,texts,autotexts = ax.pie(
    counts.values,
    colors=colors,
    startangle=90,

    wedgeprops=dict(
        width=0.42,
        edgecolor="white",
        linewidth=1
    ),

    autopct=lambda p:
        f"{p:.1f}%"
        if p>3 else "",

    pctdistance=1.1  # <-- ключевой параметр
)



    ###################################
    # center text
    ###################################

    ax.text(
        0,0,
        f"{total}\nmiRNAs",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold"
    )


    ###################################
    # legend with n
    ###################################

    leg_labels = [
        f"{lab} (n={n})"
        for lab,n in zip(
            labels,
            counts.values
        )
    ]

    ax.legend(
        wedges,
        leg_labels,
        loc="center left",
        bbox_to_anchor=(1.02,0.5),
        frameon=False
    )

    ax.set_title(
        "Pan-cancer recurrent dysregulation",
        fontsize=14,
        fontweight="bold"
    )

    plt.tight_layout()

    if save_pdf:
        plt.savefig(
            save_pdf,
            bbox_inches="tight"
        )

    plt.show()


def plot_combined_specificity_distribution(
    gtex_df: pd.DataFrame,
    gtex_col: str,
    tcga_df: pd.DataFrame,
    tcga_col: str,
    title: str = f"{SPECIFICITY_METRIC} score distribution: GTEx vs TCGA",
    save_path: str = None,
):
    # Настраиваем чистый стиль
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))

    # Палитра: приглушенный зеленый (sage/muted green) и приглушенный красный (muted red/coral)
    gtex_color = "#4c9173"  # Приглушенный зеленый
    tcga_color = "#c35151"  # Приглушенный красный

    gtex_median = gtex_df[gtex_col].median()
    tcga_median = tcga_df[tcga_col].median()

    # 1. Распределение для GTEx (Здоровые ткани)
    sns.kdeplot(
        data=gtex_df,
        x=gtex_col,
        fill=True,
        color=gtex_color,
        alpha=0.4,  # Делаем прозрачнее, чтобы области пересечения были видны
        linewidth=2,
        bw_adjust=0.5,
        label=f"GTEx (Normal) [Median: {gtex_median:.3f}]",
    )

    # 2. Распределение для TCGA (Опухоли)
    sns.kdeplot(
        data=tcga_df,
        x=tcga_col,
        fill=True,
        color=tcga_color,
        alpha=0.4,
        linewidth=2,
        bw_adjust=0.5,
        label=f"TCGA (Tumors) [Median: {tcga_median:.3f}]",
    )

    # Добавляем "ковры" (rugplot) снизу для каждого датасета с соответствующим цветом
    sns.rugplot(
        data=gtex_df, x=gtex_col, color=gtex_color, alpha=0.05, height=0.03
    )
    sns.rugplot(
        data=tcga_df, x=tcga_col, color=tcga_color, alpha=0.05, height=0.03
    )

    # Оформление осей и заголовка
    plt.title(title, fontsize=16, pad=15, fontweight="bold")
    plt.xlabel(f"{SPECIFICITY_METRIC} score", fontsize=12, labelpad=10)
    plt.ylabel("Density", fontsize=12, labelpad=10)

    plt.xlim(-0.05, 1.05)

    plt.legend(
        loc="upper left",
        fontsize=13,
        frameon=True,
        facecolor="white",
        edgecolor="gray",
    )

    sns.despine(left=True, bottom=True)
    if save_path:
        plt.savefig(
            save_path,
            bbox_inches="tight"
        )

    plt.tight_layout()
    plt.show()




def plot_specificity_scatter(
    df,
    title_suffix="",
    save_path=None,
    point_size=100,
    alpha=0.8,
):

    plot_df = df.dropna(
        subset=[f"{SPECIFICITY_METRIC}_gtex", f"{SPECIFICITY_METRIC}_tcga", "Specificity_class"]
    ).copy()


    fig, ax = plt.subplots(
        figsize=(8.2, 7.4),
        dpi=300
    )

    fig.patch.set_facecolor("white")

    # используем фиксированный порядок
    categories = [
        c for c in FINAL_CHOICES
        if c in plot_df["Specificity_class"].values
    ]

    for cat in categories:
        sub = plot_df[
            plot_df["Specificity_class"] == cat
        ]

        ax.scatter(
            sub[f"{SPECIFICITY_METRIC}_gtex"],
            sub[f"{SPECIFICITY_METRIC}_tcga"],
            s=point_size,
            alpha=alpha,
            c=[FINAL_COLORS.get(cat, "gray")],
            edgecolors="black",
            linewidths=0.4,
            label=f"{cat} (n={len(sub)})"
        )

    _specificity_lo = -SPECIFICITY_AXIS_PAD
    _specificity_hi = 1.0 + SPECIFICITY_AXIS_PAD

    # diagonal
    ax.plot(
        [_specificity_lo, _specificity_hi], [_specificity_lo, _specificity_hi],
        ls="--",
        lw=1,
        alpha=0.5
    )

    ax.set_xlabel(f"{SPECIFICITY_METRIC} GTEx", fontsize=20)
    ax.set_ylabel(f"{SPECIFICITY_METRIC} TCGA", fontsize=20)

    ax.set_title(
        f"TCGA vs GTEx {SPECIFICITY_METRIC}\n"
        + (f"\n{title_suffix}" if title_suffix else ""),
        fontsize=13
    )

    _set_specificity_scatter_axes(ax)

    ax.legend(
        bbox_to_anchor=(1.03, 1),
        loc="upper left",
        frameon=False,
        fontsize=9,
        title="Specificity class"
    )

    plt.tight_layout()

    if save_path:
        plt.savefig(
            save_path,
            bbox_inches="tight"
        )

    plt.show()


def plot_dysregulated_specificity_scatter(df, size=90, linewidths=0.4, save_path=None):
    fig, ax = plt.subplots(figsize=(6,6))

    sns.scatterplot(
        data=df,
        x=f'{SPECIFICITY_METRIC}_gtex',
        y=f'{SPECIFICITY_METRIC}_tcga',
        hue='Specificity_class',
        style='Dysregulated_class',
        edgecolors="black",
        linewidths=linewidths,
        palette=FINAL_COLORS,
        s=size,
        ax=ax
    )

    _specificity_lo = -SPECIFICITY_AXIS_PAD
    _specificity_hi = 1.0 + SPECIFICITY_AXIS_PAD

    # diagonal
    ax.plot(
        [_specificity_lo, _specificity_hi],
        [_specificity_lo, _specificity_hi],
        '--',
        color='gray',
        alpha=0.7
    )

    _set_specificity_scatter_axes(ax, equal_aspect=False)

    ax.set_xlabel(f'{SPECIFICITY_METRIC} specificity in GTEx')
    ax.set_ylabel(f'{SPECIFICITY_METRIC} specificity in TCGA')
    ax.grid(True)
    ax.set_title('Dysregulated miRNAs')
    if save_path:
        plt.savefig(save_path, bbox_inches="tight")
    plt.tight_layout()
    plt.show()



def plot_dysregulated_specificity_scatter_binary(
    df,
    title_suffix="",
    save_path=None,
    point_size=50,
    alpha=0.8,
):
    """
    Static GTEx vs TCGA specificity scatter colored by Dysregulated_class
    (same view as the last page of plot_specificity_scatter_interactive).
    """
    plot_df = df.dropna(
        subset=[
            f"{SPECIFICITY_METRIC}_gtex",
            f"{SPECIFICITY_METRIC}_tcga",
            "Dysregulated_class",
        ]
    ).copy()

    dys_order = ["UPREGULATED", "DOWNREGULATED", "NOT_SIGNIFICANT"]
    dys_labels = {
        "UPREGULATED": "Upregulated",
        "DOWNREGULATED": "Downregulated",
        "NOT_SIGNIFICANT": "Non-significant"
    }
    categories = [c for c in dys_order if c in plot_df["Dysregulated_class"].values]

    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    fig.patch.set_facecolor("white")

    for cat in categories:
        sub = plot_df[plot_df["Dysregulated_class"] == cat]
        alpha = 0.2 if cat == 'NOT_SIGNIFICANT' else alpha
        ax.scatter(
            sub[f"{SPECIFICITY_METRIC}_gtex"],
            sub[f"{SPECIFICITY_METRIC}_tcga"],
            s=point_size,
            alpha=alpha,
            c=[STN_CLASS_COLORS.get(cat, "#888888")],
            edgecolors="black",
            linewidths=0.4,
            label=f"{dys_labels.get(cat, cat)} (n={len(sub)})",
        )

    _specificity_lo = -SPECIFICITY_AXIS_PAD
    _specificity_hi = 1.0 + SPECIFICITY_AXIS_PAD

    ax.plot(
        [_specificity_lo, _specificity_hi], [_specificity_lo, _specificity_hi],
        ls="--", lw=1, alpha=0.5,
    )

    ax.set_xlabel(f"{SPECIFICITY_METRIC} GTEx", fontsize=10)
    ax.set_ylabel(f"{SPECIFICITY_METRIC} TCGA", fontsize=10)
    ax.set_title(
        f"TCGA vs GTEx {SPECIFICITY_METRIC} — Dysregulated"
        + (f"\n{title_suffix}" if title_suffix else ""),
        fontsize=10,
    )

    _set_specificity_scatter_axes(ax)

    ax.legend(
        bbox_to_anchor=(1.03, 1),
        loc="upper left",
        frameon=False,
        fontsize=6,
        title="Dysregulated class",
    )
    ax.grid(True)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, bbox_inches="tight")

    plt.show()


def plot_log2fc_specificity_scatter(
    df,
    title_suffix="",
    save_path=None,
    point_size=100,
    alpha=0.8,
    color_pct: float = 95,
    log2fc_col: str = "log2FC_tcga_gtex",
):
    """
    Static GTEx vs TCGA specificity scatter colored by log2FC (Viridis).
    Same subset as the log2FC page in plot_specificity_scatter_interactive (direction ≠ NS).
    """
    if log2fc_col not in df.columns:
        raise ValueError(f"Column {log2fc_col!r} not found in dataframe.")
    if "direction" not in df.columns:
        raise ValueError("Column 'direction' not found in dataframe.")

    plot_df = df.dropna(
        subset=[f"{SPECIFICITY_METRIC}_gtex", f"{SPECIFICITY_METRIC}_tcga", log2fc_col]
    ).copy()
    plot_df = plot_df[~plot_df["direction"].isin(["NS", "NA"])].copy()
    if plot_df.empty:
        raise ValueError("No miRNAs with direction ≠ NS and valid specificity / log2FC.")

    lim = _log2fc_color_limit(plot_df[log2fc_col], color_pct)
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)

    fig, ax = plt.subplots(figsize=(8.2, 7.4), dpi=300)
    fig.patch.set_facecolor("white")

    sc = ax.scatter(
        plot_df[f"{SPECIFICITY_METRIC}_gtex"],
        plot_df[f"{SPECIFICITY_METRIC}_tcga"],
        c=plot_df[log2fc_col],
        cmap="viridis",
        norm=norm,
        s=point_size,
        alpha=alpha,
        edgecolors="black",
        linewidths=0.4,
    )

    _specificity_lo = -SPECIFICITY_AXIS_PAD
    _specificity_hi = 1.0 + SPECIFICITY_AXIS_PAD

    ax.plot(
        [_specificity_lo, _specificity_hi], [_specificity_lo, _specificity_hi],
        ls="--", lw=1, alpha=0.5,
    )

    ax.set_xlabel(f"{SPECIFICITY_METRIC} GTEx", fontsize=20)
    ax.set_ylabel(f"{SPECIFICITY_METRIC} TCGA", fontsize=20)
    ax.set_title(
        f"TCGA vs GTEx {SPECIFICITY_METRIC} — log2FC TCGA/GTEx\n"
        f"(direction ≠ NS, n={len(plot_df)})"
        + (f"\n{title_suffix}" if title_suffix else ""),
        fontsize=13,
    )

    _set_specificity_scatter_axes(ax)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label(f"log2FC (TCGA/GTEx) [±{lim:.1f} @ p{color_pct:.0f}]", fontsize=10)

    plt.tight_layout()
    plt.grid(True)

    if save_path:
        plt.savefig(save_path, bbox_inches="tight")

    plt.show()


def plot_dysregulated_delta_vs_gtex(
    merged,
    title_suffix="",
    save_path=None,
    point_size=100,
    alpha=0.8,
    show_thresholds: bool = False,
    baseline_thr: float | None = None,
):
    """
    Static Δ specificity vs GTEx for dysregulated miRNAs only
    (same view as page 3 of plot_specificity_delta_vs_gtex).
    """
    plot_df = _prepare_specificity_delta_df(merged)
    plot_df = plot_df[
        plot_df["Dysregulated_class"].isin(["UPREGULATED", "DOWNREGULATED"])
    ].copy()
    if plot_df.empty:
        raise ValueError("No dysregulated miRNAs with valid specificity delta.")

    baseline_thr = SPECIFICITY_THR if baseline_thr is None else baseline_thr
    dys_labels = {
        "UPREGULATED": "Upregulated",
        "DOWNREGULATED": "Downregulated",
    }

    fig, ax = plt.subplots(figsize=(8.2, 6.5), dpi=300)
    fig.patch.set_facecolor("white")

    for cat in ("UPREGULATED", "DOWNREGULATED"):
        sub = plot_df[plot_df["Dysregulated_class"] == cat]
        if sub.empty:
            continue
        ax.scatter(
            sub["spec_gtex"],
            sub["spec_delta"],
            s=point_size,
            alpha=alpha,
            c=STN_CLASS_COLORS.get(cat, "#888888"),
            edgecolors="black",
            linewidths=0.4,
            label=f"{dys_labels[cat]} (n={len(sub)})",
        )

    _x_lo = -SPECIFICITY_AXIS_PAD
    _x_hi = 1.0 + SPECIFICITY_AXIS_PAD
    _y_pad = max(
    abs(plot_df["spec_delta"].min()), 
    abs(plot_df["spec_delta"].max()), 
    0.25,
    )
    _y_lo, _y_hi = -_y_pad - 0.05, _y_pad + 0.05

    ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.7)
    if show_thresholds:
        ax.axvline(baseline_thr, color="#B8860B", linewidth=1, linestyle=":", alpha=0.6)

    ax.set_xlim(_x_lo, _x_hi)
    ax.set_ylim(_y_lo, _y_hi)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel(f"{SPECIFICITY_METRIC} GTEx (normal)", fontsize=12)
    ax.set_ylabel(f"Δ {SPECIFICITY_METRIC} (TCGA − GTEx)", fontsize=12)
    ax.set_title(
        f"Δ {SPECIFICITY_METRIC} vs GTEx — dysregulated miRNAs (n={len(plot_df)})"
        + (f"\n{title_suffix}" if title_suffix else ""),
        fontsize=13,
        fontweight="bold",
    )

    ax.legend(
        bbox_to_anchor=(1.03, 1),
        loc="upper left",
        frameon=False,
        fontsize=9,
        title="Dysregulated class",
    )
    ax.grid(True, alpha=0.5)
    #sns.despine()

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, bbox_inches="tight")

    plt.show()


def _prepare_specificity_delta_df(merged: pd.DataFrame) -> pd.DataFrame:
    """Use precomputed TCGA − GTEx specificity delta; negative Δ = loss in cancer."""
    col_gtex = f"{SPECIFICITY_METRIC}_gtex"
    col_tcga = f"{SPECIFICITY_METRIC}_tcga"
    col_delta = f"Delta_{SPECIFICITY_METRIC}"
    if col_delta not in merged.columns:
        raise ValueError(f"Missing required precomputed specificity delta column: {col_delta}")
    df = merged.dropna(subset=[col_gtex, col_tcga, col_delta]).copy()
    df["spec_gtex"] = df[col_gtex].astype(float)
    df["spec_tcga"] = df[col_tcga].astype(float)
    df["spec_delta"] = df[col_delta].astype(float)
    return df.reset_index(drop=True)



def classification_plots(df_in, auc_score, save_path=None):
    # Подготовка данных
    df = df_in.copy()
    binary_df = df[df["Dysregulated_class"].isin(["DOWNREGULATED", "UPREGULATED"])].copy()
    y_true = (binary_df["Dysregulated_class"] == "UPREGULATED").astype(int)
    scores = -binary_df["Delta_Tau"]

    # Настройка стиля (минимализм)
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Цветовая схема: приглушенные, но контрастные цвета
    #colors = {"UPREGULATED": "#E67E22", "DOWNREGULATED": "#3498DB", "NOT_SIGNIFICANT": "#95A5A6"}

    # --- График 1: Плотности (KDE) ---
    for group in ["UPREGULATED", "DOWNREGULATED", "NOT_SIGNIFICANT"]:
        subset = df[df["Dysregulated_class"] == group]
        alpha=0.2 if group == 'NOT_SIGNIFICANT' else .5
        sns.kdeplot(subset["Delta_Tau"], label=group, fill=True, 
                    alpha=alpha, color=STN_CLASS_COLORS.get(group), linewidth=1, ax=ax1)

    ax1.set_title("Density of Delta_Tau Across Groups", fontsize=14, fontweight='bold', pad=15)
    ax1.set_xlabel("Delta_Tau Value", fontsize=12)
    ax1.set_ylabel("Density", fontsize=12)
    ax1.legend(frameon=False)
    sns.despine() # Убираем лишние рамки

    # --- График 2: ROC Curve ---
    fpr, tpr, _ = roc_curve(y_true, scores)
    ax2.plot(fpr, tpr, color='#2C3E50', lw=3, label=f'AUC = {auc_score:.3f}')
    ax2.plot([0, 1], [0, 1], color='#BDC3C7', linestyle='--', lw=2)
    ax2.fill_between(fpr, tpr, alpha=0.1, color='#2C3E50')
    
    ax2.set_title("ROC Analysis", fontsize=14, fontweight='bold', pad=15)
    ax2.set_xlabel("False Positive Rate", fontsize=12)
    ax2.set_ylabel("True Positive Rate", fontsize=12)
    ax2.legend(loc='lower right', frameon=False)
    sns.despine()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight")
    plt.tight_layout(pad=3.0)
    plt.show()


def _delta_sign_traces(df, hover_map, *, x_col: str, visible: bool, mk: dict) -> list:
    """Two traces per subset: loss (Δ<0) and gain (Δ≥0)."""
    traces = []
    for label, mask, color in (
        (f"Δ < 0 (loss)", df["spec_delta"] < 0, _DELTA_LOSS_COLOR),
        (f"Δ ≥ 0 (gain)", df["spec_delta"] >= 0, _DELTA_GAIN_COLOR),
    ):
        sub = df[mask]
        if sub.empty:
            continue
        traces.append(go.Scatter(
            x=sub[x_col],
            y=sub["spec_delta"],
            mode="markers",
            name=label,
            visible=visible,
            marker={**mk, "color": color},
            text=[hover_map[i] for i in sub.index],
            hoverinfo="text",
        ))
    return traces


def _dysregulated_class_traces(df, hover_map, *, x_col: str, visible: bool, mk: dict) -> list:
    """UP/DOWN traces; colors match plot_specificity_scatter_interactive."""
    traces = []
    for cat in ("UPREGULATED", "DOWNREGULATED"):
        sub = df[df["Dysregulated_class"] == cat]
        if sub.empty:
            continue
        traces.append(go.Scatter(
            x=sub[x_col],
            y=sub["spec_delta"],
            mode="markers",
            name=cat,
            visible=visible,
            marker={**mk, "color": STN_CLASS_COLORS.get(cat, "#888888")},
            text=[hover_map[i] for i in sub.index],
            hoverinfo="text",
        ))
    return traces


def plot_specificity_delta_vs_gtex(
    merged,
    save_path=None,
    show_thresholds: bool = False,
    baseline_thr: float | None = None,
    delta_thr: float | None = None,
):
    """
    Δ specificity (TCGA − GTEx) vs GTEx/TCGA baseline. Interactive Plotly.
    Views: all miRNAs vs GTEx/TCGA | dysregulated miRNAs vs GTEx/TCGA.
    """
    plot_df_all = _prepare_specificity_delta_df(merged)
    if plot_df_all.empty:
        raise ValueError("No miRNAs with valid GTEx/TCGA specificity scores.")

    plot_df_de = plot_df_all[
        plot_df_all["Dysregulated_class"].isin(["UPREGULATED", "DOWNREGULATED"])
    ].copy()

    baseline_thr = SPECIFICITY_THR if baseline_thr is None else baseline_thr
    delta_thr = delta_thr if delta_thr is not None else SPECIFICITY_AXIS_PAD * 4

    hover_all = {
        i: _mirna_hover_text(plot_df_all.loc[i], include_delta=True)
        for i in plot_df_all.index
    }
    hover_de = {
        i: _mirna_hover_text(plot_df_de.loc[i], include_delta=True)
        for i in plot_df_de.index
    }

    mk = dict(size=14, opacity=0.85, line=dict(width=0.6, color="white"))
    traces_all_gtex = _delta_sign_traces(
        plot_df_all, hover_all, x_col="spec_gtex", visible=True, mk=mk
    )
    traces_all_tcga = _delta_sign_traces(
        plot_df_all, hover_all, x_col="spec_tcga", visible=False, mk=mk
    )
    traces_de_gtex = _dysregulated_class_traces(
        plot_df_de, hover_de, x_col="spec_gtex", visible=False, mk=mk
    )
    traces_de_tcga = _dysregulated_class_traces(
        plot_df_de, hover_de, x_col="spec_tcga", visible=False, mk=mk
    )
    traces = traces_all_gtex + traces_all_tcga + traces_de_gtex + traces_de_tcga

    n_all_gtex = len(traces_all_gtex)
    n_all_tcga = len(traces_all_tcga)
    n_de_gtex = len(traces_de_gtex)
    n_de_tcga = len(traces_de_tcga)
    vis_all_gtex = (
        [True] * n_all_gtex
        + [False] * n_all_tcga
        + [False] * n_de_gtex
        + [False] * n_de_tcga
    )
    vis_all_tcga = (
        [False] * n_all_gtex
        + [True] * n_all_tcga
        + [False] * n_de_gtex
        + [False] * n_de_tcga
    )
    vis_de_gtex = (
        [False] * n_all_gtex
        + [False] * n_all_tcga
        + [True] * n_de_gtex
        + [False] * n_de_tcga
    )
    vis_de_tcga = (
        [False] * n_all_gtex
        + [False] * n_all_tcga
        + [False] * n_de_gtex
        + [True] * n_de_tcga
    )

    fig = go.Figure(data=traces)
    _x_lo = -SPECIFICITY_AXIS_PAD
    _x_hi = 1.0 + SPECIFICITY_AXIS_PAD
    _y_pad = max(
        abs(plot_df_all["spec_delta"].quantile(0.01)),
        abs(plot_df_all["spec_delta"].quantile(0.99)),
        0.25,
    )
    _y_lo, _y_hi = -_y_pad - 0.05, _y_pad + 0.05

    fig.add_shape(
        type="line", x0=_x_lo, x1=_x_hi, y0=0, y1=0,
        line=dict(color="rgba(0,0,0,0.45)", dash="dash"),
    )
    if show_thresholds:
        fig.add_shape(
            type="line", x0=baseline_thr, x1=baseline_thr, y0=_y_lo, y1=_y_hi,
            line=dict(color="rgba(184,134,11,0.6)", dash="dot"),
        )
        fig.add_shape(
            type="line", x0=_x_lo, x1=_x_hi, y0=-delta_thr, y1=-delta_thr,
            line=dict(color="rgba(100,100,100,0.35)", dash="dot"),
        )
        fig.add_shape(
            type="line", x0=_x_lo, x1=_x_hi, y0=delta_thr, y1=delta_thr,
            line=dict(color="rgba(100,100,100,0.35)", dash="dot"),
        )

    fig.update_layout(
        title=f"Δ {SPECIFICITY_METRIC} vs GTEx — all miRNAs (n={len(plot_df_all)})",
        template="plotly_white",
        height=850,
        autosize=True,
        margin=dict(l=70, r=250, t=120, b=70),
        xaxis=dict(
            title=f"{SPECIFICITY_METRIC} GTEx (normal)",
            range=[_x_lo, _x_hi],
            tickmode="array",
            tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
        ),
        yaxis=dict(
            title=f"Δ {SPECIFICITY_METRIC} (TCGA − GTEx)",
            range=[_y_lo, _y_hi],
        ),
        legend=dict(
            orientation="v",
            x=1.02, y=1,
            xanchor="left", yanchor="top",
            bordercolor="LightGrey", borderwidth=1,
        ),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0, y=1.08,
                xanchor="left", yanchor="bottom",
                buttons=[
                    dict(
                        label="All miRNA delta vs GTEx",
                        method="update",
                        args=[
                            {"visible": vis_all_gtex},
                            {
                                "title.text": (
                                f"Δ {SPECIFICITY_METRIC} vs GTEx — all miRNAs "
                                f"(n={len(plot_df_all)})"
                                ),
                                "xaxis.title.text": f"{SPECIFICITY_METRIC} GTEx (normal)",
                            },
                        ],
                    ),
                    dict(
                        label="All miRNA delta vs TCGA",
                        method="update",
                        args=[
                            {"visible": vis_all_tcga},
                            {
                                "title.text": (
                                    f"Δ {SPECIFICITY_METRIC} vs TCGA — all miRNAs "
                                    f"(n={len(plot_df_all)})"
                                ),
                                "xaxis.title.text": f"{SPECIFICITY_METRIC} TCGA (cancer)",
                            },
                        ],
                    ),
                    dict(
                        label="Dysregulated miRNAs vs GTEx",
                        method="update",
                        args=[
                            {"visible": vis_de_gtex},
                            {
                                "title.text": (
                                    f"Δ {SPECIFICITY_METRIC} vs GTEx — dysregulated miRNAs "
                                    f"(n={len(plot_df_de)})"
                                ),
                                "xaxis.title.text": f"{SPECIFICITY_METRIC} GTEx (normal)",
                            },
                        ],
                    ),
                    dict(
                        label="Dysregulated miRNAs vs TCGA",
                        method="update",
                        args=[
                            {"visible": vis_de_tcga},
                            {
                                "title.text": (
                                    f"Δ {SPECIFICITY_METRIC} vs TCGA — dysregulated miRNAs "
                                    f"(n={len(plot_df_de)})"
                                ),
                                "xaxis.title.text": f"{SPECIFICITY_METRIC} TCGA (cancer)",
                            },
                        ],
                    ),
                ],
            )
        ],
    )

    if save_path:
        fig.write_html(save_path, include_plotlyjs="cdn")
    fig.show()


def plot_specificity_delta_kde(
    merged,
    save_path=None,
    hue: str | None = "Dysregulated_class",
    title: str | None = None,
):
    """KDE of Δ specificity (TCGA − GTEx); optional split by hue column."""
    plot_df = _prepare_specificity_delta_df(merged)
    if plot_df.empty:
        raise ValueError("No miRNAs with valid GTEx/TCGA specificity scores.")

    sns.set_theme(style="whitegrid", palette="muted")
    fig, ax = plt.subplots(figsize=(9, 5))

    hue_order = None
    palette = None
    if hue and hue in plot_df.columns:
        pref = ["UPREGULATED", "DOWNREGULATED", "NOT_SIGNIFICANT"]
        seen = plot_df[hue].dropna().astype(str).unique().tolist()
        hue_order = [c for c in pref if c in seen] + sorted(c for c in seen if c not in pref)
        palette = [STN_CLASS_COLORS.get(c, "#888888") for c in hue_order]
        for cat, color in zip(hue_order, palette):
            sub = plot_df[plot_df[hue].astype(str) == cat]
            if sub.empty:
                continue
            sns.kdeplot(
                data=sub, x="spec_delta", fill=True, alpha=0.35,
                linewidth=2, bw_adjust=0.55, color=color, label=cat, ax=ax,
            )
        ax.legend(title=hue, frameon=False)
    else:
        sns.kdeplot(
            data=plot_df, x="spec_delta", fill=True, color="#2c7fb8",
            alpha=0.55, linewidth=2, bw_adjust=0.55, ax=ax,
        )

    sns.rugplot(data=plot_df, x="spec_delta", color="#555555", alpha=0.12, height=0.04, ax=ax)
    ax.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel(f"Δ {SPECIFICITY_METRIC} (TCGA − GTEx)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title(
        title or f"Distribution of Δ {SPECIFICITY_METRIC} (n={len(plot_df)})",
        fontsize=14, fontweight="bold", pad=12,
    )
    sns.despine(left=True, bottom=True)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.show()


def plot_specificity_delta_waterfall(
    merged,
    save_path=None,
    color_by: str = "sign",
    label_extremes: int = 0,
    figsize=(14, 6),
):
    """
    Waterfall: miRNAs sorted by Δ specificity (loss → gain).
    color_by: 'sign' | 'Specificity_class' | column name in merged.
    """
    plot_df = _prepare_specificity_delta_df(merged).sort_values("spec_delta").reset_index(drop=True)
    if plot_df.empty:
        raise ValueError("No miRNAs with valid GTEx/TCGA specificity scores.")

    n = len(plot_df)
    x = np.arange(n)

    if color_by == "sign":
        colors = np.where(plot_df["spec_delta"] < 0, "#DC143C", "#9467BD")
    elif color_by in plot_df.columns:
        cmap = {**FINAL_COLORS, **STN_CLASS_COLORS}
        colors = plot_df[color_by].astype(str).map(lambda c: cmap.get(c, "#888888"))
    else:
        colors = "#2c7fb8"

    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(x, plot_df["spec_delta"], color=colors, width=1.0, edgecolor="none", alpha=0.9)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="-", alpha=0.7)

    if label_extremes > 0:
        k = min(label_extremes, n // 2)
        label_idx = list(range(k)) + list(range(n - k, n))
        for i in label_idx:
            row = plot_df.iloc[i]
            ax.text(
                i, row["spec_delta"],
                row["mature_name"],
                ha="center",
                va="bottom" if row["spec_delta"] >= 0 else "top",
                fontsize=7,
                rotation=90,
            )

    ax.set_xlabel(f"miRNAs sorted by Δ {SPECIFICITY_METRIC} (loss ← → gain)", fontsize=11)
    ax.set_ylabel(f"Δ {SPECIFICITY_METRIC} (TCGA − GTEx)", fontsize=11)
    ax.set_title(f"Specificity change waterfall (n={n})", fontsize=13, fontweight="bold")
    ax.set_xlim(-0.5, n - 0.5)
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    sns.despine()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.show()


def plot_volcano(
    res_proj: pd.DataFrame,
    ax,
    title: str,
    *,
    padj_max: float = PADJ_THR,
    abs_lfc_min: float = LFC_THR,
    label_top_padj: int = 0,
    label_top_abs_lfc: int = 0,
) -> None:
    """Volcano: up=firebrick, down=navy, faint gray between thresholds; labels optional."""
    r = res_proj.copy()
    r["padj"] = pd.to_numeric(r["padj"], errors="coerce")
    r["log2FoldChange"] = pd.to_numeric(r["log2FoldChange"], errors="coerce")
    padj = np.clip(r["padj"].to_numpy(dtype=float), 1e-300, None)
    r["neglog10padj"] = -np.log10(padj)
    lfc = r["log2FoldChange"]

    ok = r["padj"].notna() & lfc.notna() & r["neglog10padj"].notna()
    hi = ok & (r["padj"] < padj_max)
    up = hi & (lfc > abs_lfc_min)
    dn = hi & (lfc < -abs_lfc_min)
    between = ok & ~up & ~dn

    ax.set_facecolor("white")
    ax.figure.patch.set_facecolor("white")

    if between.any():
        ax.scatter(
            r.loc[between, "log2FoldChange"],
            r.loc[between, "neglog10padj"],
            s=24,
            alpha=0.22,
            c="#888888",
            edgecolors="#aaaaaa",
            linewidths=0.35,
            rasterized=True,
            zorder=1,
        )
    if dn.any():
        ax.scatter(
            r.loc[dn, "log2FoldChange"],
            r.loc[dn, "neglog10padj"],
            s=26,
            alpha=0.92,
            c="navy",
            edgecolors="black",
            linewidths=0.55,
            rasterized=True,
            zorder=3,
        )
    if up.any():
        ax.scatter(
            r.loc[up, "log2FoldChange"],
            r.loc[up, "neglog10padj"],
            s=26,
            alpha=0.92,
            c="firebrick",
            edgecolors="black",
            linewidths=0.55,
            rasterized=True,
            zorder=3,
        )

    y_thr = -np.log10(padj_max)
    ax.axhline(y_thr, color="#444444", linestyle="--", linewidth=1.0, zorder=0)
    ax.axvline(abs_lfc_min, color="#444444", linestyle=":", linewidth=0.95, zorder=0)
    ax.axvline(-abs_lfc_min, color="#444444", linestyle=":", linewidth=0.95, zorder=0)

    ax.text(
        0.02,
        0.98,
        f"guides: FDR={padj_max}, |LFC|={abs_lfc_min}",
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#cccccc", alpha=0.92),
    )

    if "mature_name" in r.columns and (label_top_padj > 0 or label_top_abs_lfc > 0):
        r_ok = r.loc[ok].dropna(subset=["mature_name"])
        if not r_ok.empty:
            chunks = []
            r_up = r_ok[r_ok["log2FoldChange"] >= 0]
            r_dn = r_ok[r_ok["log2FoldChange"] < 0]

            if label_top_padj > 0:
                if not r_up.empty:
                    chunks.append(r_up.nsmallest(label_top_padj, "padj"))
                if not r_dn.empty:
                    chunks.append(r_dn.nsmallest(label_top_padj, "padj"))

            if label_top_abs_lfc > 0:
                if not r_up.empty:
                    chunks.append(
                        r_up.assign(abs_lfc=r_up["log2FoldChange"].abs()).nlargest(
                            label_top_abs_lfc, "abs_lfc"
                        )
                    )
                if not r_dn.empty:
                    chunks.append(
                        r_dn.assign(abs_lfc=r_dn["log2FoldChange"].abs()).nlargest(
                            label_top_abs_lfc, "abs_lfc"
                        )
                    )

            to_label = pd.concat(chunks).drop_duplicates(subset=["mature_name"])

            xs, ys, cols = [], [], []
            for _, row in to_label.iterrows():
                x = float(row["log2FoldChange"])
                y = float(row["neglog10padj"])
                pj = row["padj"]
                xs.append(x)
                ys.append(y)
                if pd.notna(pj) and pj < padj_max and x > abs_lfc_min:
                    cols.append("firebrick")
                elif pd.notna(pj) and pj < padj_max and x < -abs_lfc_min:
                    cols.append("navy")
                else:
                    cols.append("#555555")
            ax.scatter(
                xs,
                ys,
                s=28,
                c=cols,
                edgecolors="black",
                linewidths=0.6,
                alpha=0.95,
                zorder=6,
                rasterized=True,
            )

            for _, row in to_label.iterrows():
                x, y = float(row["log2FoldChange"]), float(row["neglog10padj"])
                if x < 0:
                    xytext, ha = (-10, 6), "right"
                else:
                    xytext, ha = (10, 6), "left"
                ax.annotate(
                    str(row["mature_name"]),
                    xy=(x, y),
                    xytext=xytext,
                    textcoords="offset points",
                    fontsize=7,
                    ha=ha,
                    arrowprops=dict(arrowstyle="-", color="0.25", lw=0.65, shrinkA=0, shrinkB=2),
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", alpha=0.95, linewidth=0.6),
                    zorder=7,
                )

    ax.set_xlabel("log2FC (tumor vs normal)", fontsize=10)
    ax.set_ylabel(r"$-\log_{10}$ padj", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="600")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)





def plot_pan_cancer_clustermap(
    fc_mat,
    row_colors,
    row_cluster=False,
    col_cluster=True,
    save_pdf=None
):

    vmax = max(
        3,
        np.nanpercentile(
            np.abs(fc_mat.values),
            98
        )
    )

    # fixed scale centered so
    # -1..1 becomes beige plateau
    vmin = -vmax


    cg = sns.clustermap(
        fc_mat,

        cmap=custom_cmap,

        vmin=vmin,
        vmax=vmax,

        row_cluster=row_cluster,
        col_cluster=col_cluster,

        metric="correlation",
        method="average",

        row_colors=row_colors,

        figsize=(12,14),

        linewidths=.25,
        linecolor="black",

        xticklabels=True,
        yticklabels=True,

        dendrogram_ratio=(.1,.12),

        cbar_pos=(
            .02,.82,.03,.12
        ),

        cbar_kws={
            "label":"log2FC tumor vs normal"
        }
    )


    ##################################
    # labels
    ##################################

    cg.ax_heatmap.set_xlabel(
        "TCGA projects",
        fontsize=12
    )

    cg.ax_heatmap.set_ylabel(
        "Recurrent miRNAs",
        fontsize=12
    )

    cg.fig.suptitle(
        "Pan-cancer recurrent miRNA dysregulation",
        y=1.02,
        fontsize=15,
        fontweight="bold"
    )


    ##################################
    # legend
    ##################################

    for lab,col in {
        "UP":"firebrick",
        "DOWN":"navy"
    }.items():

        cg.ax_col_dendrogram.bar(
            0,0,
            color=col,
            label=lab,
            linewidth=0
        )

    cg.ax_col_dendrogram.legend(
        frameon=False,
        ncol=2,
        title="Binomial class"
    )


    if save_pdf:
        cg.savefig(
            save_pdf,
            bbox_inches="tight"
        )

    plt.show()


_EXPRESSION_HEATMAP_CMAPS = {"viridis", "plasma", "inferno", "rocket"}


def _expression_heatmap_matrix(
    gtex_cpm_aggregate: pd.DataFrame,
    tcga_cpm_aggregate: pd.DataFrame,
    mirnas: list[str] | tuple[str, ...] | pd.Series,
) -> pd.DataFrame:
    mirnas = pd.Index(pd.Series(mirnas, dtype="object").dropna().astype(str)).drop_duplicates()
    if mirnas.empty:
        raise ValueError("mirnas must contain at least one miRNA name.")

    present = [
        mir for mir in mirnas
        if mir in gtex_cpm_aggregate.columns and mir in tcga_cpm_aggregate.columns
    ]
    missing_any = [mir for mir in mirnas if mir not in present]
    if not present:
        raise ValueError("None of the requested miRNAs were found in both GTEx and TCGA columns.")
    if missing_any:
        warnings.warn(
            "Some miRNAs were skipped because they were not found in both datasets: "
            + ", ".join(missing_any),
            stacklevel=2,
        )

    parts = []
    for label, df in (("GTEx", gtex_cpm_aggregate), ("TCGA", tcga_cpm_aggregate)):
        mat = df.loc[:, present].astype(float).T
        mat.columns = [f"{label}: {col}" for col in mat.columns.astype(str)]
        parts.append(mat)

    return pd.concat(parts, axis=1)


def _apply_heatmap_z_score(mat: pd.DataFrame, z_score: str | bool | None) -> tuple[pd.DataFrame, str]:
    if z_score in (None, False, "none", "None"):
        return mat, "log2(CPM + 1)"

    if z_score is True:
        z_score = "row"
    if z_score not in {"row", "column"}:
        raise ValueError("z_score must be one of None, 'row', 'column', or True.")

    if z_score == "row":
        mean = mat.mean(axis=1)
        std = mat.std(axis=1, ddof=0).replace(0, np.nan)
        return mat.sub(mean, axis=0).div(std, axis=0).fillna(0), "row z-score"

    mean = mat.mean(axis=0)
    std = mat.std(axis=0, ddof=0).replace(0, np.nan)
    return mat.sub(mean, axis=1).div(std, axis=1).fillna(0), "column z-score"


def _sort_heatmap_axis(mat: pd.DataFrame, axis: int, mode: str | bool | None) -> pd.DataFrame:
    if mode in (None, False, "none", "None", "cluster"):
        return mat

    if mode is True:
        mode = "mean"
    if mode not in {"mean", "median", "max", "name"}:
        raise ValueError("sort_rows/sort_cols must be False, True, 'mean', 'median', 'max', 'name', or 'cluster'.")

    if mode == "name":
        labels = sorted(mat.index if axis == 0 else mat.columns)
    elif axis == 0:
        labels = getattr(mat, mode)(axis=1).sort_values(ascending=False).index
    else:
        labels = getattr(mat, mode)(axis=0).sort_values(ascending=False).index

    return mat.loc[labels, :] if axis == 0 else mat.loc[:, labels]


def plot_mirna_expression_heatmap(
    gtex_cpm_aggregate: pd.DataFrame,
    tcga_cpm_aggregate: pd.DataFrame,
    mirnas: list[str] | tuple[str, ...] | pd.Series,
    *,
    z_score: str | bool | None = "row",
    sort_rows: str | bool | None = False,
    sort_cols: str | bool | None = False,
    cmap: str = "plasma",
    figsize: tuple[float, float] | None = None,
    linewidths: float = 0.5,
    linecolor: str = "black",
    title: str | None = None,
    save_path: str | None = None,
    show: bool = True,
):
    """
    Plot miRNA expression heatmap from aggregated log2(CPM+1) GTEx and TCGA tables.

    Rows are miRNAs, columns are GTEx tissues and TCGA projects.
    z_score='row' highlights tissue/project specificity patterns; use z_score=None
    to show absolute log2(CPM+1) expression.

    sort_rows/sort_cols: False, True/'mean', 'median', 'max', 'name', or 'cluster'.
    cmap: 'plasma', 'viridis', 'inferno', or 'rocket'.
    """
    if cmap not in _EXPRESSION_HEATMAP_CMAPS:
        raise ValueError(f"cmap must be one of: {', '.join(sorted(_EXPRESSION_HEATMAP_CMAPS))}.")

    mat = _expression_heatmap_matrix(gtex_cpm_aggregate, tcga_cpm_aggregate, mirnas)
    mat, cbar_label = _apply_heatmap_z_score(mat, z_score)
    mat = _sort_heatmap_axis(mat, axis=0, mode=sort_rows)
    mat = _sort_heatmap_axis(mat, axis=1, mode=sort_cols)

    row_cluster = sort_rows == "cluster"
    col_cluster = sort_cols == "cluster"
    if figsize is None:
        figsize = (
            max(10, min(28, 0.32 * mat.shape[1] + 4)),
            max(5, min(24, 0.35 * mat.shape[0] + 3)),
        )

    sns.set_theme(style="white")
    cg = sns.clustermap(
        mat,
        cmap=cmap,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        linewidths=linewidths,
        linecolor=linecolor,
        figsize=figsize,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": cbar_label},
        dendrogram_ratio=(0.08, 0.12), 
    )

    cg.ax_heatmap.set_xlabel("GTEx tissues / TCGA projects", fontsize=11)
    cg.ax_heatmap.set_ylabel("miRNA", fontsize=11)
    cg.ax_heatmap.tick_params(axis="x", labelrotation=90, labelsize=8)
    cg.ax_heatmap.tick_params(axis="y", labelsize=9)
    cg.fig.suptitle(
        title or f"miRNA expression heatmap ({cbar_label})",
        y=1.02,
        fontsize=14,
        fontweight="bold",
    )

    if save_path:
        cg.savefig(save_path, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return cg



def plot_boxplot_gtex_tcga(df_1, 
                           df_2, 
                           mir_name=None,
                           labels=('GTEx', 'TCGA'), 
                           ylabel='Expression Level', 
                           figsize=(8, 4),
                           save_path: str = None):

    if mir_name not in df_1.columns or mir_name not in df_2.columns:
        raise ValueError(f"miRNA '{mir_name}' not found in the dataframe columns.")
    
    x = df_1[mir_name].to_numpy()
    y = df_2[mir_name].to_numpy()
    
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    
    # 2. Считаем статистику (Mann-Whitney U-test)
    if len(x) > 0 and len(y) > 0:
        _, pval = mannwhitneyu(x, y, alternative='two-sided')
        if pval < 0.001:
            p_text = "p < 0.001"
        elif pval < 0.05:
            p_text = f"p = {pval:.3f}"
        else:
            p_text = f"p = {pval:.3f} (ns)"
    else:
        p_text = "p = NaN"

    # 3. Собираем длинный датафрейм для Seaborn
    df_plot = pd.DataFrame({
        'Expression': np.concatenate([x, y]),
        'Group': [labels[0]] * len(x) + [labels[1]] * len(y)
    })
    
    # 4. Настройка темы графика
    sns.set_theme(style="whitegrid", rc={"axes.facecolor": "#f9f9f9"})
    fig, ax = plt.subplots(figsize=figsize)
    
    # Контрастная палитра для двух групп (SteelBlue и IndianRed)
    palette = {labels[0]: "#4682B4", labels[1]: "#CD5C5C"}
    
    # Рисуем горизонтальные боксплоты
    sns.boxplot(
        data=df_plot,
        x='Expression',
        y='Group',
        hue='Group',
        palette=palette,
        legend=False,
        fliersize=0,
        width=0.4,
        linewidth=1.2,
        ax=ax
    )
    
    # Накладываем точки с тонкой черной обводкой
    sns.stripplot(
        data=df_plot,
        x='Expression',
        y='Group',
        hue='Group',
        palette=palette,
        legend=False,
        size=4,
        alpha=0.5,
        jitter=0.15,
        linewidth=0.6,
        edgecolor='black',
        ax=ax
    )
    
    # 5. Кастомизация подписей
    # Поскольку график перевернут, шкала экспрессии (ylabel) становится осью X
    ax.set_xlabel(ylabel, fontsize=11, labelpad=10)
    ax.set_ylabel("", fontsize=11) # Названия групп на оси Y говорят сами за себя
    
    # Формируем красивый заголовок с p-value
    full_title = f"{mir_name} ({p_text})" if mir_name else p_text
    ax.set_title(full_title, fontsize=13, fontweight='bold', pad=15)
    
    # Настройка сетки: только вертикальные линии
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, linestyle='--', alpha=0.7)
    
    # Убираем рамки
    sns.despine(left=True, bottom=True)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.tight_layout()
    return fig, ax


def plot_mirna_boxplots_gtex_tcga(
    df_1: pd.DataFrame,
    df_2: pd.DataFrame,
    mirnas: list[str] | tuple[str, ...] | pd.Series,
    labels: tuple[str, str] = ("GTEx", "TCGA"),
    colors: tuple[str, str] | list[str] | dict[str, str] = ("#4682B4", "#CD5C5C"),
    ylabel: str = "log2(CPM+1)",
    figsize: tuple[float, float] | None = None,
    max_quantile: float = 1.0,
    show_points: bool = True,
    point_size: float = 3,
    point_alpha: float = 0.45,
    save_path: str | None = None,
    show: bool = True,
):
    """
    Horizontal panel of vertical GTEx vs TCGA expression boxplots for several miRNAs.

    df_1 and df_2 are expected to be aggregated log2(CPM+1) tables with miRNAs
    in columns and tissues/projects in rows. No statistical test is calculated.
    """
    mirnas = pd.Index(pd.Series(mirnas, dtype="object").dropna().astype(str)).drop_duplicates()
    if mirnas.empty:
        raise ValueError("mirnas must contain at least one miRNA name.")

    if isinstance(colors, dict):
        palette = {label: colors[label] for label in labels}
    else:
        if len(colors) != 2:
            raise ValueError("colors must contain exactly two colors or be a dict keyed by labels.")
        palette = dict(zip(labels, colors))

    present = [mir for mir in mirnas if mir in df_1.columns and mir in df_2.columns]
    missing = [mir for mir in mirnas if mir not in present]
    if not present:
        raise ValueError("None of the requested miRNAs were found in both dataframes.")
    if missing:
        warnings.warn(
            "Some miRNAs were skipped because they were not found in both dataframes: "
            + ", ".join(missing),
            stacklevel=2,
        )

    records = []
    for mir in present:
        for label, df in zip(labels, (df_1, df_2)):
            values = pd.to_numeric(df[mir], errors="coerce").dropna()
            records.extend(
                {"miRNA": mir, "Group": label, "Expression": value}
                for value in values
            )

    plot_df = pd.DataFrame(records)
    if plot_df.empty:
        raise ValueError("No non-NA expression values found for requested miRNAs.")

    n_mirnas = len(present)
    if figsize is None:
        figsize = (max(6, 2.35 * n_mirnas), 3)

    sns.set_theme(
        style="whitegrid",
        rc={
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "grid.color": "#e6e6e6",
            "grid.linewidth": 0.8,
        },
    )
    fig, axes = plt.subplots(1, n_mirnas, figsize=figsize, sharey=True)
    axes = np.atleast_1d(axes)

    finite_expr = plot_df["Expression"].replace([np.inf, -np.inf], np.nan).dropna()
    y_min = min(0, finite_expr.min())
    y_max = finite_expr.quantile(max_quantile)
    if np.isclose(y_min, y_max):
        y_max = y_min + 1

    for ax, mir in zip(axes, present):
        sub = plot_df[plot_df["miRNA"] == mir]
        sns.boxplot(
            data=sub,
            x="Group",
            y="Expression",
            hue="Group",
            order=list(labels),
            hue_order=list(labels),
            palette=palette,
            legend=False,
            fliersize=0,
            width=0.58,
            linewidth=1.25,
            saturation=0.95,
            boxprops={"edgecolor": "#222222", "linewidth": 1.25, "alpha": 0.9},
            whiskerprops={"color": "#222222", "linewidth": 1.15},
            capprops={"color": "#222222", "linewidth": 1.15},
            medianprops={"color": "black", "linewidth": 1.9},
            ax=ax,
        )
        if show_points:
            sns.stripplot(
                data=sub,
                x="Group",
                y="Expression",
                hue="Group",
                order=list(labels),
                hue_order=list(labels),
                palette=palette,
                legend=False,
                size=point_size,
                alpha=point_alpha,
                jitter=0.15,
                linewidth=0.5,
                edgecolor="black",
                zorder=3,
                ax=ax,
            )

        ax.set_title(mir, fontsize=12, fontweight="bold", pad=8)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_ylim(y_min, y_max)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True, linestyle="--", alpha=0.75)
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=10)
        sns.despine(ax=ax, left=False, bottom=True)

    axes[0].set_ylabel(ylabel, fontsize=11, labelpad=10)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes



def plot_expression_boxplots(df: pd.DataFrame, 
                            gene_name: str, 
                            tissue_col: str = 'tissue', 
                            figsize=(8, 12),
                            max_quantile: float = 1,
                            save_path: str = None):
    # Проверяем наличие гена в данных
    if gene_name not in df.columns:
        raise ValueError(f"Gen '{gene_name}' not found in the dataframe columns.")
    
    # Вытаскиваем нужные данные  
    if tissue_col == 'index':
        plot_data = df[[gene_name]].copy()
        plot_data['tissue_clean'] = df.index
    else:
        if tissue_col not in df.columns:
            raise ValueError(f"Column with tissues '{tissue_col}' not found in the dataframe.")
        plot_data = df[[tissue_col, gene_name]].copy()
        plot_data.rename(columns={tissue_col: 'tissue_clean'}, inplace=True)
        
    # Очищаем данные от NaN в экспрессии гена
    plot_data = plot_data.dropna(subset=[gene_name])
    
    # cортируем ткани по убыванию медианы (высокая экспрессия будет вверху само собой)
    tissue_order = (plot_data.groupby('tissue_clean')[gene_name]
                    .median()
                    .sort_values(ascending=False)
                    .index)
    
    # стиль
    sns.set_theme(style="whitegrid", rc={"axes.facecolor": "white"})

    
    ########
    tissue_order = (plot_data.groupby('tissue_clean')[gene_name]
                    .median()
                    .sort_values(ascending=False)
                    .index)
    palette_colors = sns.color_palette("viridis", n_colors=len(tissue_order))[::-1]
    color_dict = dict(zip(tissue_order, palette_colors))
    sns.set_theme(style="whitegrid", rc={"axes.facecolor": "white"})
    fig, ax = plt.subplots(figsize=figsize)
    
    # Рисуем горизонтальные боксплоты: x — экспрессия, y — ткань
    sns.boxplot(
        data=plot_data, 
        x=gene_name, 
        y='tissue_clean', 
        order=tissue_order,
        hue_order=tissue_order,
        palette=color_dict,
        fliersize=0,
        width=0.6,
        linewidth=1.2,
        ax=ax
    )
    
    # Накладываем индивидуальные точки (Stripplot)
    sns.stripplot(
        data=plot_data, 
        x=gene_name, 
        y='tissue_clean', 
        order=tissue_order,
        color="#DB7093", 
        size=1.5, 
        alpha=0.7, 
        jitter=0.15,
        dodge=False,
        ax=ax
    )
    limit = plot_data[gene_name].quantile(max_quantile)
    ax.set_xlim(0, limit)
    
    # Кастомизация подписей и эстетика
    ax.set_title(f"{gene_name}", fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel("Expression Level", fontsize=11, labelpad=10)
    ax.set_ylabel("Tissue / Cell Type", fontsize=11, labelpad=10)
    
    # Настройки шрифтов для осей (названия тканей теперь строго горизонтальные на оси Y)
    ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='x', labelsize=10)
    
    # Переключаем сетку на вертикальную (по оси X)
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, linestyle='--', alpha=0.7)
    
    # Убираем лишние рамки
    sns.despine(left=True, bottom=True)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.tight_layout()
    plt.show()