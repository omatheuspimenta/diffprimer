"""
HTML Report Generator for diffprimer.

Produces a standalone, interactive HTML report with Plotly.js charts,
a searchable DataTable, genome browser visualization, and full parameter log.
"""

import os
import re
import json
import html as html_mod
import base64
from pathlib import Path

import pandas as pd
from diffprimer import __version__
from diffprimer.logs import diffprimerLog

logger = diffprimerLog()

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _encode_logo_base64() -> str:
    """Read the logo PNG and return a data-URI string.  Falls back to empty."""
    logo_paths = [
        Path(__file__).resolve().parents[2] / "docs" / "logo.png",          # dev layout
        Path(__file__).resolve().parents[3] / "docs" / "logo.png",          # installed
    ]
    for p in logo_paths:
        if p.exists():
            raw = p.read_bytes()
            b64 = base64.b64encode(raw).decode("ascii")
            return f"data:image/png;base64,{b64}"
    return ""


def _parse_header(header_col):
    """
    Parses the header column formatted as:
    OriginalHeader_region_N_Start_End
    Returns a DataFrame with Sequence_ID, Region_ID, Start, End.
    """
    pattern = re.compile(r"(.*)_(region_\d+)_(\d+)_(\d+)$")
    parsed = []
    for h in header_col:
        match = pattern.match(str(h))
        if match:
            parsed.append(match.groups())
        else:
            parsed.append((h, "Unknown", "Unknown", "Unknown"))
    return pd.DataFrame(parsed, columns=["Sequence_ID", "Region_ID", "Start", "End"])

# ──────────────────────────────────────────────────────────────────────────────
# Color palette — single source of truth
# ──────────────────────────────────────────────────────────────────────────────
PALETTE = {
    "primary":        "#2563eb",
    "primary_light":  "#3b82f6",
    "primary_dark":   "#1d4ed8",
    "success":        "#059669",
    "success_light":  "#34d399",
    "danger":         "#dc2626",
    "danger_light":   "#f87171",
    "warning":        "#d97706",
    "warning_light":  "#fbbf24",
    "purple":         "#7c3aed",
    "purple_light":   "#a78bfa",
    "teal":           "#0d9488",
    "teal_light":     "#5eead4",
    "slate":          "#64748b",
    "bg":             "#f8fafc",
    "surface":        "#ffffff",
    "border":         "#e2e8f0",
    "text":           "#1e293b",
    "text_muted":     "#64748b",
}

CHART_COLORS = [
    PALETTE["primary"], PALETTE["success"], PALETTE["purple"],
    PALETTE["warning"], PALETTE["teal"], PALETTE["danger"],
    PALETTE["primary_light"], PALETTE["success_light"],
]

# ──────────────────────────────────────────────────────────────────────────────
# Tooltip descriptions for every table column
# ──────────────────────────────────────────────────────────────────────────────
COLUMN_TOOLTIPS = {
    "Sequence_ID":         "Identifier of the reference contig where the exclusive region was found.",
    "Region_ID":           "Sequential identifier for the exclusive region within this contig.",
    "Start":               "Genomic start coordinate (bp) of the exclusive region in the reference.",
    "End":                 "Genomic end coordinate (bp) of the exclusive region in the reference.",
    "Forw_Seq":            "Nucleotide sequence of the forward (left) primer (5' → 3').",
    "Rev_Seq":             "Nucleotide sequence of the reverse (right) primer (5' → 3').",
    "Left_Size":           "Length (nt) of the forward primer.",
    "Right_Size":          "Length (nt) of the reverse primer.",
    "Amplicon_Seq":        "Full nucleotide sequence of the predicted amplicon.",
    "Amplicon_Size":       "Length (bp) of the predicted PCR amplicon.",
    "Region_Length":       "Total length (bp) of the exclusive region.",
    "Annotation":          "Overlapping genomic feature annotations from the GFF3 file.",
    "Specificity_Tag":     "Classification of primer specificity: Specific (low global similarity or positional mismatches) vs Non-Specific.",
    "Most_Similar_Target": "Off-target sequence with the highest similarity to the amplicon region.",
    "Max_Sim%":            "Percentage similarity of the amplicon region to the most similar off-target sequence.",
    "Max Sim%":            "Percentage similarity of the amplicon region to the most similar off-target sequence.",
}

# ──────────────────────────────────────────────────────────────────────────────
# Main generator
# ──────────────────────────────────────────────────────────────────────────────

def generate_html_report(
    results_dir: str,
    run_params: dict | None = None,
    primer3_args: dict | None = None,
    contig_info: list[dict] | None = None,
    version: str | None = None,
):
    """Generate a standalone interactive HTML report.

    Parameters
    ----------
    results_dir : str
        Directory containing ``primer_design_results.csv``.
    run_params : dict, optional
        Pipeline parameters (reference_file, k, etc.).
    primer3_args : dict, optional
        Primer3 global arguments used in the run.
    contig_info : list[dict], optional
        ``[{header, length, regions: [{start, end}]}]`` for genome viz.
    version : str, optional
        Software version string.
    """

    csv_file = os.path.join(results_dir, "primer_design_results.csv")
    if not os.path.exists(csv_file):
        logger.error(f"Cannot generate report: {csv_file} does not exist.")
        return

    logger.info("Generating HTML report...")

    # ── Load data ─────────────────────────────────────────────────────────
    try:
        df = pd.read_csv(csv_file, sep=";")
    except Exception as e:
        logger.error(f"Failed to read CSV file: {e}")
        return

    # Parse header into Sequence_ID / Region_ID / Start / End
    parsed_header_df = _parse_header(df["Header"])

    # Ensure expected columns exist
    expected_cols = [
        "Header", "Forw_Seq", "Rev_Seq", "Amplicon_Seq", "Amplicon_Size",
        "Region_Length", "Annotation", "Sequence_Region", "Specificity_Tag",
        "Most_Similar_Target", "Max_Similarity_pct",
    ]
    for col in expected_cols:
        if col not in df.columns:
            df[col] = "NA"

    table_df = pd.concat([parsed_header_df, df[expected_cols]], axis=1)

    # Compute primer sizes
    table_df["Left_Size"] = df["Forw_Seq"].apply(lambda x: len(str(x)) if str(x) != "NA" else "NA")
    table_df["Right_Size"] = df["Rev_Seq"].apply(lambda x: len(str(x)) if str(x) != "NA" else "NA")

    # ── Pre-compute data for charts ───────────────────────────────────────
    amp_sizes = df["Amplicon_Size"].replace("NA", pd.NA).dropna().astype(float).tolist()

    spec_tags = df["Specificity_Tag"].value_counts().to_dict()
    spec_labels = list(spec_tags.keys())
    spec_values = list(spec_tags.values())

    left_tm  = df.get("Left_TM",  pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    right_tm = df.get("Right_TM", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    left_gc  = df.get("Left_GC",  pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    right_gc = df.get("Right_GC", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()

    left_hairpin = df.get("Left_HAIRPIN_TH",  pd.Series()).replace(["NA", "-"], pd.NA).dropna().astype(float).tolist()
    right_hairpin= df.get("Right_HAIRPIN_TH", pd.Series()).replace(["NA", "-"], pd.NA).dropna().astype(float).tolist()
    left_self    = df.get("Left_SELF_ANY_TH", pd.Series()).replace(["NA", "-"], pd.NA).dropna().astype(float).tolist()
    right_self   = df.get("Right_SELF_ANY_TH",pd.Series()).replace(["NA", "-"], pd.NA).dropna().astype(float).tolist()

    max_sim_pct = df.get("Max_Similarity_pct", pd.Series()).replace(["NA", "-"], pd.NA).dropna().astype(float).tolist()

    # ── Logo ──────────────────────────────────────────────────────────────
    logo_data_uri = _encode_logo_base64()

    # ── Version ───────────────────────────────────────────────────────────
    sw_version = version or __version__

    # ── Contig info for genome viz ────────────────────────────────────────
    # Enrich contig_info with primer positions and annotations for the Genome Map
    if contig_info:
        for contig in contig_info:
            for r in contig.get("regions", []):
                # Match region by Sequence_ID and Start
                c_head = contig.get("header")
                r_start = str(r.get("start"))
                match = table_df[(table_df["Sequence_ID"] == c_head) & (table_df["Start"].astype(str) == r_start)]
                if not match.empty:
                    row = match.iloc[0]
                    # Parse gene ID for genome map tooltip
                    ann_raw = str(row.get("Annotation", ""))
                    gene_info = ann_raw
                    for part in ann_raw.replace("|", ";").split(";"):
                        if part.startswith("gene=") or part.startswith("gene_id=") or part.startswith("Name=") or part.startswith("ID="):
                            gene_info = part.split("=")[1]
                            break
                    r["annotation"] = gene_info
                    
                    seq_reg = str(row.get("Sequence_Region", ""))
                    forw_seq = str(row.get("Forw_Seq", ""))
                    if forw_seq and forw_seq not in ("NA", "NaN") and forw_seq in seq_reg:
                        idx = seq_reg.index(forw_seq)
                        r["primer_left_start"] = idx
                        r["primer_left_len"] = len(forw_seq)
                        
                        amp_s = str(row.get("Amplicon_Size", "0"))
                        if amp_s.replace(".", "").isdigit():
                            amp_size = int(float(amp_s))
                            rev_seq = str(row.get("Rev_Seq", ""))
                            if amp_size > 0 and rev_seq and rev_seq not in ("NA", "NaN"):
                                r["primer_right_len"] = len(rev_seq)
                                r["primer_right_start"] = idx + amp_size - len(rev_seq)

    contig_json = json.dumps(contig_info or [])

    # ── Specificity counts ────────────────────────────────────────────────
    n_designed   = len(amp_sizes)
    n_specific   = spec_tags.get("Specific_HighSim", 0) + spec_tags.get("Specific_LowGlobalSim", 0) + spec_tags.get("Specific_PositionalMismatches", 0)
    n_nonspecific = spec_tags.get("NonSpecific_Amplification", 0)
    n_not_checked = spec_tags.get("Not_Checked", 0)

    # ── Build the table columns & rows ────────────────────────────────────
    # Columns in display order (Header removed, Left/Right Size added)
    display_columns = [
        "Sequence_ID", "Region_ID", "Start", "End",
        "Forw_Seq", "Rev_Seq", "Left_Size", "Right_Size",
        "Amplicon_Seq", "Amplicon_Size", "Region_Length",
        "Annotation", "Specificity_Tag", "Most_Similar_Target", "Max_Similarity_pct",
    ]
    display_labels = {
        "Sequence_ID": "Sequence ID",
        "Region_ID": "Region ID",
        "Start": "Start",
        "End": "End",
        "Forw_Seq": "Forward Primer",
        "Rev_Seq": "Reverse Primer",
        "Left_Size": "Fwd Size",
        "Right_Size": "Rev Size",
        "Amplicon_Seq": "Amplicon Seq",
        "Amplicon_Size": "Amplicon Size",
        "Region_Length": "Region Length",
        "Annotation": "Annotation",
        "Specificity_Tag": "Specificity",
        "Most_Similar_Target": "Most Similar Target",
        "Max_Similarity_pct": "Max Sim%",
    }

    # Build header HTML with tooltips
    table_header_html = ""
    for col in display_columns:
        label = display_labels.get(col, col)
        tip = COLUMN_TOOLTIPS.get(col, COLUMN_TOOLTIPS.get(label, ""))
        tooltip_attr = f' title="{tip}"' if tip else ""
        help_icon = f' <span class="th-help" {tooltip_attr}>?</span>' if tip else ""
        table_header_html += f"<th>{label}{help_icon}</th>\n"

    # Build rows HTML
    def _copy_cell(full_value, display_value, truncate=False):
        fv_safe = html_mod.escape(str(full_value))
        dv_safe = html_mod.escape(str(display_value))
        if truncate:
            display_html = f'<span class="seq-cell" title="{fv_safe}">{dv_safe}</span>'
        else:
            display_html = f'<span class="copy-text">{dv_safe}</span>'
        fv_escaped = str(full_value).replace("\\", "\\\\").replace("'", "\\'")
        icon = '<svg fill="currentColor" viewBox="0 0 20 20"><path d="M8 3a1 1 0 011-1h2a1 1 0 110 2H9a1 1 0 01-1-1z"></path><path d="M6 3a2 2 0 00-2 2v11a2 2 0 002 2h8a2 2 0 002-2V5a2 2 0 00-2-2 3 3 0 01-3 3H9a3 3 0 01-3-3z"></path></svg>'
        return f'<div class="flex-copy-cell">{display_html} <button class="copy-btn" onclick="copyText(this, \'{fv_escaped}\')" title="Copy">{icon}</button></div>'

    table_body_html = ""
    for _, row in table_df.iterrows():
        tag = str(row.get("Specificity_Tag", "Not_Checked"))
        badge_cls = "notchecked"
        if "Specific" in tag and "Non" not in tag:
            badge_cls = "specific"
        elif "NonSpecific" in tag:
            badge_cls = "nonspecific"

        tag_html = f'<span class="badge {badge_cls}">{tag}</span>'

        cells = ""
        for col in display_columns:
            val = row.get(col, "")
            if val == "NA":
                val = "-"
            
            if col == "Annotation":
                # Preprocess output for HTML text only
                val_str = str(val)
                if val_str not in ("-", "no annotation", "no annotation (file not loaded)"):
                    gene_info = val_str
                    for part in val_str.replace("|", ";").split(";"):
                        if part.startswith("gene=") or part.startswith("gene_id=") or part.startswith("Name=") or part.startswith("ID="):
                            gene_info = part.split("=")[1]
                            break
                    val = gene_info

            if col == "Specificity_Tag":
                cells += f"<td>{tag_html}</td>"
            elif col in ("Forw_Seq", "Rev_Seq", "Amplicon_Seq", "Annotation", "Most_Similar_Target"):
                cells += f"<td>{_copy_cell(val, val, truncate=True)}</td>"
            else:
                cells += f"<td>{_copy_cell(val, val)}</td>"

        table_body_html += f"<tr>{cells}</tr>\n"

    # ── Run parameters table HTML ─────────────────────────────────────────
    params_html = ""
    if run_params:
        params_html += '<table class="params-table"><thead><tr><th>Parameter</th><th>Value</th></tr></thead><tbody>'
        for k_param, v_param in run_params.items():
            params_html += f"<tr><td>{k_param}</td><td>{v_param}</td></tr>"
        params_html += "</tbody></table>"

    p3_html = ""
    if primer3_args:
        p3_html += '<table class="params-table"><thead><tr><th>Parameter</th><th>Value</th></tr></thead><tbody>'
        for k_param, v_param in primer3_args.items():
            p3_html += f"<tr><td>{k_param}</td><td>{v_param}</td></tr>"
        p3_html += "</tbody></table>"

    # ── Assemble full HTML ────────────────────────────────────────────────
    logo_img_tag = f'<img src="{logo_data_uri}" alt="diffprimer logo" class="header-logo">' if logo_data_uri else ""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>diffprimer Report</title>
    <meta name="description" content="Interactive PCR primer design report generated by diffprimer v{sw_version}">
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
    <link rel="stylesheet" href="https://cdn.datatables.net/2.0.3/css/dataTables.dataTables.css" />
    <style>
/* ── CSS Reset & Variables ─────────────────────────────────────────── */
:root {{
    --primary: {PALETTE["primary"]};
    --primary-light: {PALETTE["primary_light"]};
    --primary-dark: {PALETTE["primary_dark"]};
    --success: {PALETTE["success"]};
    --success-light: {PALETTE["success_light"]};
    --danger: {PALETTE["danger"]};
    --danger-light: {PALETTE["danger_light"]};
    --warning: {PALETTE["warning"]};
    --warning-light: {PALETTE["warning_light"]};
    --purple: {PALETTE["purple"]};
    --teal: {PALETTE["teal"]};
    --bg: {PALETTE["bg"]};
    --surface: {PALETTE["surface"]};
    --border: {PALETTE["border"]};
    --text: {PALETTE["text"]};
    --text-muted: {PALETTE["text_muted"]};
    --shadow-sm: 0 1px 2px rgba(0,0,0,0.05);
    --shadow: 0 4px 6px -1px rgba(0,0,0,0.07), 0 2px 4px -2px rgba(0,0,0,0.05);
    --shadow-lg: 0 10px 15px -3px rgba(0,0,0,0.08), 0 4px 6px -4px rgba(0,0,0,0.05);
    --radius: 12px;
    --radius-sm: 8px;
}}

*, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

body {{
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.6;
    -webkit-font-smoothing: antialiased;
}}

/* ── Navigation ────────────────────────────────────────────────────── */
.top-nav {{
    position: sticky;
    top: 0;
    z-index: 100;
    background: rgba(255,255,255,0.85);
    backdrop-filter: blur(12px);
    border-bottom: 1px solid var(--border);
    padding: 0 5%;
    transition: box-shadow 0.3s;
}}
.top-nav.scrolled {{
    box-shadow: var(--shadow);
}}
.nav-inner {{
    max-width: 1400px;
    margin: 0 auto;
    display: flex;
    align-items: center;
    gap: 2rem;
    height: 52px;
    overflow-x: auto;
}}
.nav-brand {{
    font-weight: 700;
    font-size: 1rem;
    color: var(--primary);
    white-space: nowrap;
    display: flex;
    align-items: center;
    gap: 0.5rem;
}}
.nav-brand img {{ height: 28px; }}
.nav-links {{
    display: flex;
    gap: 0.25rem;
    list-style: none;
}}
.nav-links a {{
    text-decoration: none;
    color: var(--text-muted);
    font-size: 0.82rem;
    font-weight: 500;
    padding: 0.4rem 0.75rem;
    border-radius: 6px;
    transition: all 0.2s;
    white-space: nowrap;
}}
.nav-links a:hover,
.nav-links a.active {{
    color: var(--primary);
    background: rgba(37,99,235,0.08);
}}

/* ── Header ────────────────────────────────────────────────────────── */
.hero {{
    background: linear-gradient(135deg, {PALETTE["primary_dark"]} 0%, {PALETTE["primary"]} 40%, {PALETTE["primary_light"]} 100%);
    color: white;
    padding: 3rem 5%;
    position: relative;
    overflow: hidden;
}}
.hero::before {{
    content: '';
    position: absolute;
    top: -50%;
    right: -20%;
    width: 600px;
    height: 600px;
    background: radial-gradient(circle, rgba(255,255,255,0.08) 0%, transparent 70%);
    border-radius: 50%;
}}
.hero-inner {{
    max-width: 1400px;
    margin: 0 auto;
    position: relative;
    z-index: 2;
    display: flex;
    align-items: center;
    gap: 2rem;
}}
.header-logo {{
    height: 72px;
    filter: brightness(0) invert(1);
    opacity: 0.95;
}}
.hero h1 {{
    font-size: 2.2rem;
    font-weight: 700;
    letter-spacing: -0.5px;
    margin-bottom: 0.25rem;
}}
.hero-subtitle {{
    font-size: 1.05rem;
    opacity: 0.85;
    font-weight: 400;
}}
.hero-version {{
    display: inline-block;
    margin-top: 0.5rem;
    background: rgba(255,255,255,0.15);
    padding: 0.2rem 0.75rem;
    border-radius: 20px;
    font-size: 0.8rem;
    font-weight: 500;
}}

/* ── Container ─────────────────────────────────────────────────────── */
.container {{
    max-width: 1400px;
    margin: 0 auto;
    padding: 2rem 5%;
}}

/* ── Section Titles ────────────────────────────────────────────────── */
.section-title {{
    font-size: 1.4rem;
    font-weight: 700;
    color: var(--text);
    margin-bottom: 1.25rem;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--border);
    display: flex;
    align-items: center;
    gap: 0.5rem;
}}
.section-title .icon {{
    width: 24px; height: 24px;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    background: var(--primary);
    color: white;
    border-radius: 6px;
    font-size: 0.75rem;
    font-weight: 700;
}}

/* ── Stat Cards ────────────────────────────────────────────────────── */
.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 1.25rem;
    margin-bottom: 2.5rem;
}}

.stat-card {{
    background: var(--surface);
    border-radius: var(--radius);
    padding: 1.5rem;
    border: 1px solid var(--border);
    box-shadow: var(--shadow-sm);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
    animation: fadeUp 0.5s ease both;
}}
.stat-card:nth-child(1) {{ animation-delay: 0.05s; }}
.stat-card:nth-child(2) {{ animation-delay: 0.1s; }}
.stat-card:nth-child(3) {{ animation-delay: 0.15s; }}
.stat-card:nth-child(4) {{ animation-delay: 0.2s; }}

.stat-card:hover {{
    transform: translateY(-3px);
    box-shadow: var(--shadow-lg);
}}
.stat-card h3 {{
    margin: 0 0 0.5rem 0;
    font-size: 0.78rem;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.8px;
    font-weight: 600;
}}
.stat-card .value {{
    font-size: 2.4rem;
    font-weight: 700;
    line-height: 1;
}}
.stat-card .value.blue   {{ color: var(--primary); }}
.stat-card .value.green  {{ color: var(--success); }}
.stat-card .value.red    {{ color: var(--danger); }}
.stat-card .value.amber  {{ color: var(--warning); }}

@keyframes fadeUp {{
    from {{ opacity: 0; transform: translateY(20px); }}
    to   {{ opacity: 1; transform: translateY(0); }}
}}

/* ── Charts Grid ───────────────────────────────────────────────────── */
.charts-grid {{
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 1.5rem;
    margin-bottom: 2.5rem;
}}
@media (max-width: 900px) {{
    .charts-grid {{ grid-template-columns: 1fr; }}
}}

.chart-card {{
    background: var(--surface);
    border-radius: var(--radius);
    border: 1px solid var(--border);
    box-shadow: var(--shadow-sm);
    overflow: hidden;
    animation: fadeUp 0.5s ease both;
}}
.chart-card.full-width {{
    grid-column: 1 / -1;
}}
.chart-header {{
    padding: 1.25rem 1.5rem 0;
}}
.chart-header h2 {{
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--text);
    margin: 0;
}}
.chart-header p {{
    font-size: 0.82rem;
    color: var(--text-muted);
    margin: 0.35rem 0 0;
    line-height: 1.5;
}}
.chart-body {{
    padding: 0.75rem 1rem 1rem;
}}

/* ── Table ─────────────────────────────────────────────────────────── */
.table-section {{
    background: var(--surface);
    border-radius: var(--radius);
    border: 1px solid var(--border);
    box-shadow: var(--shadow-sm);
    padding: 1.5rem 2rem 2rem;
    margin-bottom: 2.5rem;
    overflow-x: auto;
    animation: fadeUp 0.5s ease both;
}}
table.dataTable {{
    border-collapse: collapse !important;
    width: 100% !important;
    color: var(--text);
}}
table.dataTable th, table.dataTable td {{
    font-family: 'Inter', sans-serif;
    font-size: 0.82rem;
    padding: 10px 12px;
    border-bottom: 1px solid var(--border);
}}
table.dataTable thead th {{
    background-color: rgba(0,0,0,0.02);
    color: var(--text-muted);
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.4px;
    font-size: 0.72rem;
    border-bottom: 2px solid var(--border);
    white-space: nowrap;
}}
table.dataTable tbody tr {{
    background-color: transparent !important;
    transition: background-color 0.15s;
}}
table.dataTable tbody tr:hover {{
    background-color: rgba(37,99,235,0.03) !important;
}}
.dt-top {{
    display: flex;
    align-items: center;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: 0.75rem;
    margin-bottom: 1rem;
}}
.dt-search {{
    display: flex;
    align-items: center;
    gap: 0.75rem;
    flex-wrap: wrap;
}}
.dt-search input {{
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 8px 12px;
    background: var(--surface);
    color: var(--text);
    font-family: 'Inter', sans-serif;
    font-size: 0.82rem;
}}
.custom-dt-filter {{
    display: inline-flex;
    align-items: center;
    gap: 6px;
    font-size: 0.82rem;
    color: var(--text-muted);
    font-weight: 500;
}}
.custom-dt-filter select {{
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 8px 12px;
    background: var(--surface);
    color: var(--text);
    font-family: 'Inter', sans-serif;
    font-size: 0.82rem;
    cursor: pointer;
    outline: none;
    transition: border-color 0.2s;
}}
.custom-dt-filter select:focus {{
    border-color: var(--primary-light);
    box-shadow: 0 0 0 3px rgba(37,99,235,0.1);
}}
.th-help {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 16px; height: 16px;
    border-radius: 50%;
    background: var(--primary);
    color: white;
    font-size: 0.6rem;
    font-weight: 700;
    margin-left: 4px;
    cursor: help;
    vertical-align: middle;
    opacity: 0.7;
    transition: opacity 0.2s;
}}
.th-help:hover {{ opacity: 1; }}

.seq-cell {{
    font-family: 'Inter', sans-serif;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    max-width: 140px;
    display: inline-block;
    vertical-align: middle;
    font-size: 0.82rem;
}}
.flex-copy-cell {{
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 4px;
}}
.copy-text {{
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}}
.copy-btn {{
    background: none;
    border: none;
    cursor: pointer;
    color: var(--text-muted);
    margin-left: 4px;
    padding: 2px;
    opacity: 0.4;
    transition: opacity 0.2s, color 0.2s;
    vertical-align: middle;
    display: inline-flex;
    align-items: center;
}}
.copy-btn:hover {{ opacity: 1; color: var(--primary); }}
.copy-btn svg {{ width: 13px; height: 13px; }}

.badge {{
    display: inline-block;
    padding: 3px 8px;
    border-radius: 20px;
    font-size: 0.7rem;
    font-weight: 600;
    letter-spacing: 0.3px;
}}
.badge.specific    {{ background: #dcfce7; color: #166534; }}
.badge.nonspecific  {{ background: #fee2e2; color: #991b1b; }}
.badge.notchecked  {{ background: #f1f5f9; color: #475569; }}

/* ── Parameters Footer ─────────────────────────────────────────────── */
.params-section {{
    background: var(--surface);
    border-radius: var(--radius);
    border: 1px solid var(--border);
    box-shadow: var(--shadow-sm);
    margin-bottom: 2.5rem;
    overflow: hidden;
    animation: fadeUp 0.5s ease both;
}}
.params-toggle {{
    width: 100%;
    padding: 1rem 1.5rem;
    background: none;
    border: none;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: space-between;
    font-family: 'Inter', sans-serif;
    font-size: 1rem;
    font-weight: 600;
    color: var(--text);
    transition: background 0.2s;
}}
.params-toggle:hover {{ background: rgba(0,0,0,0.02); }}
.params-toggle .chevron {{
    transition: transform 0.3s;
    font-size: 1.2rem;
    color: var(--text-muted);
}}
.params-toggle.open .chevron {{
    transform: rotate(180deg);
}}
.params-body {{
    max-height: 0;
    overflow: hidden;
    transition: max-height 0.4s ease;
}}
.params-body.open {{
    max-height: 2000px;
}}
.params-body-inner {{
    padding: 0 1.5rem 1.5rem;
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
    gap: 1.5rem;
}}
.params-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.82rem;
}}
.params-table th {{
    text-align: left;
    padding: 8px 12px;
    background: rgba(0,0,0,0.02);
    border-bottom: 2px solid var(--border);
    color: var(--text-muted);
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.4px;
    font-size: 0.72rem;
}}
.params-table td {{
    padding: 7px 12px;
    border-bottom: 1px solid var(--border);
    color: var(--text);
    word-break: break-all;
}}
.params-table tbody tr:hover {{
    background: rgba(37,99,235,0.03);
}}

/* ── Footer ────────────────────────────────────────────────────────── */
.footer {{
    text-align: center;
    padding: 2rem;
    color: var(--text-muted);
    font-size: 0.8rem;
    border-top: 1px solid var(--border);
}}
.footer a {{
    color: var(--primary);
    text-decoration: none;
}}

/* ── Utility ───────────────────────────────────────────────────────── */
.mb-2 {{ margin-bottom: 2rem; }}
.mt-1 {{ margin-top: 1rem; }}

    </style>
</head>
<body>

<!-- ─── Top Navigation ─────────────────────────────────────────────── -->
<nav class="top-nav" id="topNav">
    <div class="nav-inner">
        <div class="nav-brand">
            {f'<img src="{logo_data_uri}" alt="logo">' if logo_data_uri else ''}
            diffprimer
        </div>
        <ul class="nav-links">
            <li><a href="#summary">Summary</a></li>
            <li><a href="#charts">Charts</a></li>
            <li><a href="#genome-map">Genome Map</a></li>
            <li><a href="#data-table">Data Table</a></li>
            <li><a href="#parameters">Parameters</a></li>
        </ul>
    </div>
</nav>

<!-- ─── Hero Header ────────────────────────────────────────────────── -->
<section class="hero">
    <div class="hero-inner">
        {logo_img_tag}
        <div>
            <h1>diffprimer Report</h1>
            <div class="hero-subtitle">Interactive Analysis of Target-Specific PCR Primers</div>
            <span class="hero-version">v{sw_version}</span>
        </div>
    </div>
</section>

<div class="container">

    <!-- ─── Summary Stats ──────────────────────────────────────────── -->
    <section id="summary" class="mb-2">
        <h2 class="section-title"><span class="icon">Σ</span> Summary</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <h3>Designed Pairs</h3>
                <div class="value blue">{n_designed}</div>
            </div>
            <div class="stat-card">
                <h3>Specific Pairs</h3>
                <div class="value green">{n_specific}</div>
            </div>
            <div class="stat-card">
                <h3>Non-Specific Pairs</h3>
                <div class="value red">{n_nonspecific}</div>
            </div>
            <div class="stat-card">
                <h3>Not Checked</h3>
                <div class="value amber">{n_not_checked}</div>
            </div>
        </div>
    </section>

    <!-- ─── Charts ─────────────────────────────────────────────────── -->
    <section id="charts" class="mb-2">
        <h2 class="section-title"><span class="icon">◍</span> Visualizations</h2>
        <div class="charts-grid">
            <!-- Amplicon Size -->
            <div class="chart-card">
                <div class="chart-header">
                    <h2>Amplicon Size Distribution</h2>
                    <p>Distribution of predicted PCR amplicon lengths (bp).</p>
                </div>
                <div class="chart-body"><div id="plot-amplicon"></div></div>
            </div>

            <!-- Specificity Donut -->
            <div class="chart-card">
                <div class="chart-header">
                    <h2>Specificity Breakdown</h2>
                    <p>Proportion of primers classified by specificity analysis.<br/><br/>
                    <strong style="font-size:0.75rem;">Specific_LowGlobalSim</strong>: No off-target hit above similarity threshold.<br/>
                    <strong style="font-size:0.75rem;">Specific_PositionalMismatches</strong>: Similarity found, but 3' mismatches prevent amplification.<br/>
                    <strong style="font-size:0.75rem;">NonSpecific_Amplification</strong>: Will amplify an off-target sequence.</p>
                </div>
                <div class="chart-body"><div id="plot-specificity"></div></div>
            </div>

            <!-- Max Homology -->
            <div class="chart-card">
                <div class="chart-header">
                    <h2>Off-Target Homology</h2>
                    <p>Distribution of maximum percentage similarity between each amplicon region and the closest off-target sequence. Lower off-target homology percentages indicate that the selected amplicon is highly unique to the target region, reducing the chance of cross-reacting amplification in other parts of the genome.</p>
                </div>
                <div class="chart-body"><div id="plot-homology"></div></div>
            </div>

            <!-- Dimerization -->
            <div class="chart-card">
                <div class="chart-header">
                    <h2>Dimerization &amp; Hairpin Risks</h2>
                    <p>Thermodynamic propensity scores computed by Primer3 for secondary structures.
                       <strong>Hairpin</strong>: tendency of a single primer to fold onto itself.
                       <strong>Self-Any</strong>: tendency of two copies of the same primer to dimerize.
                       Lower scores indicate lower risk. Values above 45°C may need attention.</p>
                </div>
                <div class="chart-body"><div id="plot-dimer"></div></div>
            </div>

            <!-- Tm vs GC -->
            <div class="chart-card full-width">
                <div class="chart-header">
                    <h2>Primer Parameter Space (Tm vs GC)</h2>
                    <p>Each point represents a primer. Tight clustering indicates consistent thermodynamic properties across the designed set.</p>
                </div>
                <div class="chart-body"><div id="plot-tmgc"></div></div>
            </div>
        </div>
    </section>

    <!-- ─── Genome Map ─────────────────────────────────────────────── -->
    <section id="genome-map" class="mb-2">
        <div class="params-section">
            <button class="params-toggle open" onclick="toggleParams(this, 'genomeVizBody')">
                <span>⊞ Genome Map</span>
                <span class="chevron">▾</span>
            </button>
            <div class="params-body open" id="genomeVizBody" style="padding: 1rem 1.5rem;">
                <div style="margin-bottom: 1rem; display: flex; justify-content: flex-end;">
                    <label style="font-size: 0.85rem; color: var(--text); display: flex; align-items: center; gap: 0.5rem; cursor: pointer;">
                        <input type="checkbox" id="toggleAnnotations" checked onchange="updateGenomeAnnotations(this.checked)" />
                        Show Genome Annotations
                    </label>
                </div>
                <div id="genome-viz-container"></div>
            </div>
        </div>
    </section>

    <!-- ─── Data Table ─────────────────────────────────────────────── -->
    <section id="data-table" class="mb-2">
        <h2 class="section-title"><span class="icon">☰</span> Detailed Primer Data</h2>
        <div class="table-section">
            <table id="primerTable" class="display" style="width:100%">
                <thead><tr>
                    {table_header_html}
                </tr></thead>
                <tbody>
                    {table_body_html}
                </tbody>
            </table>
        </div>
    </section>

    <!-- ─── Parameters ─────────────────────────────────────────────── -->
    <section id="parameters" class="mb-2">
        <div class="params-section">
            <button class="params-toggle" onclick="toggleParams(this, 'paramsBody')">
                <span>⚙ Run Parameters &amp; Software Information</span>
                <span class="chevron">▾</span>
            </button>
            <div class="params-body" id="paramsBody">
                <div class="params-body-inner">
                    <div>
                        <h3 style="font-size:0.95rem; margin-bottom:0.75rem; color:var(--text);">Pipeline Parameters</h3>
                        {params_html if params_html else '<p style="color:var(--text-muted); font-size:0.85rem;">No parameter data available.</p>'}
                    </div>
                    <div>
                        <h3 style="font-size:0.95rem; margin-bottom:0.75rem; color:var(--text);">Primer3 Configuration</h3>
                        {p3_html if p3_html else '<p style="color:var(--text-muted); font-size:0.85rem;">No Primer3 configuration data available.</p>'}
                    </div>
                </div>
                <div style="padding: 0 1.5rem 1rem; font-size:0.8rem; color:var(--text-muted);">
                    diffprimer v{sw_version} &mdash; Report generated automatically at the end of the analysis pipeline.
                </div>
            </div>
        </div>
    </section>

</div>

<!-- ─── Footer ─────────────────────────────────────────────────────── -->
<footer class="footer">
    Generated by <a href="https://github.com/omatheuspimenta/diffprimer" target="_blank">diffprimer</a> v{sw_version}
</footer>

<!-- ─── Scripts ────────────────────────────────────────────────────── -->
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/2.0.3/js/dataTables.js"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>

<script>
// ── Utility ──────────────────────────────────────────────────────────
function copyText(button, text) {{
    navigator.clipboard.writeText(text).then(() => {{
        const orig = button.innerHTML;
        button.innerHTML = '<svg fill="currentColor" viewBox="0 0 20 20" style="color:#059669;"><path fill-rule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clip-rule="evenodd"></path></svg>';
        setTimeout(() => {{ button.innerHTML = orig; }}, 1500);
    }});
}}

function toggleParams(btn, bodyId) {{
    btn.classList.toggle('open');
    if (!bodyId) {{ bodyId = 'paramsBody'; }}
    document.getElementById(bodyId).classList.toggle('open');
}}

// Map Plotly SVG markers to circles
function roundPlotlyLegends() {{
    document.querySelectorAll('g.legend path').forEach(path => {{
        const d = path.getAttribute('d');
        // Replace typical square path with a small circle path
        if (d && d.includes('Z') && d.length > 5 && d.length < 50) {{
            path.setAttribute('d', 'M 0,-4 A 4,4 0 1,1 0,4 A 4,4 0 1,1 0,-4 Z');
        }}
    }});
}}
setInterval(roundPlotlyLegends, 1000);

// ── Navigation scroll spy ────────────────────────────────────────────
const nav = document.getElementById('topNav');
window.addEventListener('scroll', () => {{
    nav.classList.toggle('scrolled', window.scrollY > 10);
}});

// Smooth scroll
document.querySelectorAll('.nav-links a').forEach(a => {{
    a.addEventListener('click', e => {{
        e.preventDefault();
        const target = document.querySelector(a.getAttribute('href'));
        if (target) target.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
    }});
}});

// Active section highlight
const sections = document.querySelectorAll('section[id]');
const navLinks = document.querySelectorAll('.nav-links a');
const observer = new IntersectionObserver(entries => {{
    entries.forEach(entry => {{
        if (entry.isIntersecting) {{
            navLinks.forEach(l => l.classList.remove('active'));
            const active = document.querySelector(`.nav-links a[href="#${{entry.target.id}}"]`);
            if (active) active.classList.add('active');
        }}
    }});
}}, {{ threshold: 0.3 }});
sections.forEach(s => observer.observe(s));

// ── DataTable ────────────────────────────────────────────────────────
$(document).ready(function () {{
    var table = $('#primerTable').DataTable({{
        pageLength: 25,
        scrollX: true,
        language: {{ search: "Filter:" }},
        dom: '<"dt-top"lf>rt<"dt-bottom"ip>',
        initComplete: function () {{
            var api = this.api();
            // Find the Specificity column
            api.columns().every(function () {{
                var column = this;
                var headerText = $(column.header()).text().trim().replace('?', '').trim();
                
                if (headerText === 'Specificity') {{
                    // Build filter widget and place it in the toolbar (next to search box)
                    var wrapper = $('<div class="custom-dt-filter"></div>');
                    var label = $('<span>Specificity:</span>');
                    var select = $('<select><option value="">All</option></select>')
                        .on('change', function () {{
                            var val = $.fn.dataTable.util.escapeRegex($(this).val());
                            column.search(val ? '^' + val + '$' : '', true, false).draw();
                        }});
                    // Populate unique values
                    column.data().unique().sort().each(function (d) {{
                        var textVal = $(d).text();
                        if (textVal) select.append('<option value="' + textVal + '">' + textVal + '</option>');
                    }});
                    wrapper.append(label).append(select);
                    // Append to the search/filter toolbar area
                    var searchContainer = $(api.table().container()).find('.dt-search');
                    if (searchContainer.length) {{
                        searchContainer.append(wrapper);
                    }} else {{
                        $(api.table().container()).find('.dt-top').append(wrapper);
                    }}
                }}
            }});
        }}
    }});
}});

// ── Plotly shared config ─────────────────────────────────────────────
const COLORS = {json.dumps(CHART_COLORS)};
const layoutBase = {{
    font: {{ family: 'Inter, sans-serif', size: 12, color: '{PALETTE["text"]}' }},
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    margin: {{ t: 30, r: 30, l: 55, b: 55 }},
    hoverlabel: {{ font: {{ family: 'Inter' }} }}
}};
const plotConfig = {{ responsive: true, displayModeBar: false }};

// ── Plot 1: Amplicon Size ────────────────────────────────────────────
Plotly.newPlot('plot-amplicon', [
    {{
        x: {json.dumps(amp_sizes)},
        type: 'histogram',
        marker: {{
            color: 'rgba(37,99,235,0.35)',
            line: {{ color: '{PALETTE["primary"]}', width: 1.5 }}
        }},
        name: 'Distribution',
        showlegend: false
    }},
    {{
        x: {json.dumps(amp_sizes)},
        type: 'box',
        marker: {{ color: '{PALETTE["primary"]}' }},
        line: {{ color: '{PALETTE["primary_dark"]}' }},
        fillcolor: 'rgba(37,99,235,0.15)',
        name: 'Box',
        boxmean: true,
        yaxis: 'y2',
        showlegend: false
    }}
], {{
    ...layoutBase,
    xaxis: {{ title: 'Amplicon Size (bp)', gridcolor: 'rgba(0,0,0,0.06)' }},
    yaxis: {{ title: 'Count', gridcolor: 'rgba(0,0,0,0.06)', domain: [0, 0.72] }},
    yaxis2: {{ domain: [0.78, 1], showticklabels: false, showgrid: false, zeroline: false }},
    showlegend: false,
    margin: {{ t: 20, r: 30, l: 50, b: 55 }},
    bargap: 0.05
}}, plotConfig);

// ── Plot 2: Specificity Donut ────────────────────────────────────────
Plotly.newPlot('plot-specificity', [{{
    labels: {json.dumps(spec_labels)},
    values: {json.dumps(spec_values)},
    type: 'pie',
    hole: 0.55,
    marker: {{ colors: COLORS, line: {{ color: '#fff', width: 2 }} }},
    textinfo: 'percent',
    textposition: 'outside',
    textfont: {{ size: 18 }},
    hoverinfo: 'label+value+percent',
    pull: 0.02
}}], {{
    ...layoutBase,
    showlegend: true,
    legend: {{ orientation: 'h', y: -0.15, x: 0.5, xanchor: 'center', font: {{ size: 12 }}, itemsizing: 'constant', itemwidth: 50}},
    margin: {{ t: 20, r: 40, l: 40, b: 60 }}
}}, plotConfig).then(function() {{
    setTimeout(function() {{
        var svg = document.querySelector('#plot-specificity svg');
        if (!svg) return;

        document.querySelectorAll('#plot-specificity .legendpie')
            .forEach(function(el) {{
                var fill = el.style.fill || el.getAttribute('fill') || 'gray';

                // Pega posição real na tela
                var rect = el.getBoundingClientRect();
                var svgRect = svg.getBoundingClientRect();

                // Converte para coordenadas do SVG
                var pt = svg.createSVGPoint();
                pt.x = rect.left + rect.width / 2;
                pt.y = rect.top + rect.height / 2;
                var svgP = pt.matrixTransform(svg.getScreenCTM().inverse());

                // Cria círculo no SVG raiz (fora do clipPath)
                var circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                circle.setAttribute('cx', svgP.x);
                circle.setAttribute('cy', svgP.y);
                circle.setAttribute('r', '8');
                circle.setAttribute('fill', fill);

                svg.appendChild(circle);
                el.style.display = 'none';
            }});
    }}, 400);
}});

// ── Plot 3: Off-Target Homology ──────────────────────────────────────
Plotly.newPlot('plot-homology', [{{
    x: {json.dumps(max_sim_pct)},
    type: 'histogram',
    marker: {{
        color: 'rgba(124,58,237,0.6)',
        line: {{ color: '{PALETTE["purple"]}', width: 1.5 }}
    }},
    opacity: 0.85
}}], {{
    ...layoutBase,
    xaxis: {{ title: 'Similarity to Off-Target (%)', gridcolor: 'rgba(0,0,0,0.06)' }},
    yaxis: {{ title: 'Count', gridcolor: 'rgba(0,0,0,0.06)' }}
}}, plotConfig);

// ── Plot 4: Dimerization & Hairpin ───────────────────────────────────
(function() {{
    const traces = [
        {{
            x: {json.dumps(left_hairpin)},
            y: {json.dumps(left_hairpin)}.map(() => 'Fwd Hairpin'),
            name: 'Fwd Hairpin',
            type: 'box',
            orientation: 'h',
            marker: {{ color: '{PALETTE["primary"]}' }},
            line: {{ width: 1.5 }},
            fillcolor: 'rgba(37,99,235,0.2)',
            boxmean: true
        }},
        {{
            x: {json.dumps(right_hairpin)},
            y: {json.dumps(right_hairpin)}.map(() => 'Rev Hairpin'),
            name: 'Rev Hairpin',
            type: 'box',
            orientation: 'h',
            marker: {{ color: '{PALETTE["teal"]}' }},
            line: {{ width: 1.5 }},
            fillcolor: 'rgba(13,148,136,0.2)',
            boxmean: true
        }},
        {{
            x: {json.dumps(left_self)},
            y: {json.dumps(left_self)}.map(() => 'Fwd Self-Any'),
            name: 'Fwd Self-Any',
            type: 'box',
            orientation: 'h',
            marker: {{ color: '{PALETTE["warning"]}' }},
            line: {{ width: 1.5 }},
            fillcolor: 'rgba(217,119,6,0.2)',
            boxmean: true
        }},
        {{
            x: {json.dumps(right_self)},
            y: {json.dumps(right_self)}.map(() => 'Rev Self-Any'),
            name: 'Rev Self-Any',
            type: 'box',
            orientation: 'h',
            marker: {{ color: '{PALETTE["purple"]}' }},
            line: {{ width: 1.5 }},
            fillcolor: 'rgba(124,58,237,0.2)',
            boxmean: true
        }}
    ];
    Plotly.newPlot('plot-dimer', traces, {{
        ...layoutBase,
        xaxis: {{ title: 'ΔG Propensity Score (°C)', gridcolor: 'rgba(0,0,0,0.06)' }},
        showlegend: false,
        margin: {{ t: 20, r: 30, b: 55, l: 100 }}
    }}, plotConfig);
}})();

// ── Plot 5: Tm vs GC ────────────────────────────────────────────────
Plotly.newPlot('plot-tmgc', [
    {{
        x: {json.dumps(left_gc)},
        y: {json.dumps(left_tm)},
        mode: 'markers',
        type: 'scatter',
        name: 'Forward',
        marker: {{ size: 7, color: '{PALETTE["primary"]}', opacity: 0.7, line: {{ color: 'white', width: 0.5 }} }}
    }},
    {{
        x: {json.dumps(right_gc)},
        y: {json.dumps(right_tm)},
        mode: 'markers',
        type: 'scatter',
        name: 'Reverse',
        marker: {{ size: 7, color: '{PALETTE["warning"]}', opacity: 0.7, line: {{ color: 'white', width: 0.5 }} }}
    }}
], {{
    ...layoutBase,
    xaxis: {{ title: 'GC Content (%)', gridcolor: 'rgba(0,0,0,0.06)' }},
    yaxis: {{ title: 'Melting Temperature (°C)', gridcolor: 'rgba(0,0,0,0.06)' }},
    legend: {{ orientation: 'h', y: -0.18, x: 0.5, xanchor: 'center' }}
}}, plotConfig);

// ── Genome Browser Visualization ─────────────────────────────────────
window.genomePlots = [];
window.updateGenomeAnnotations = function(show) {{
    window.genomePlots.forEach(plotId => {{
        const el = document.getElementById(plotId);
        if (el && el.layout && el.layout.annotations) {{
            const newAnns = el.layout.annotations.map(a => {{
                if (a.name === 'gene_annotation') {{
                    a.visible = show;
                }}
                return a;
            }});
            Plotly.relayout(plotId, {{ annotations: newAnns }});
        }}
    }});
}};

(function() {{
    const contigData = {contig_json};
    const container = document.getElementById('genome-viz-container');

    if (!contigData || contigData.length === 0) {{
        container.innerHTML = '<div class="chart-card" style="padding:2rem;text-align:center;color:var(--text-muted);">No contig data available for genome visualization.</div>';
        return;
    }}

    // Do not filter out any contigs to show full information for all contigs
    const contigsWithRegions = contigData;
    
    // Build one chart card per contig
    contigsWithRegions.forEach((contig, cIdx) => {{
        const card = document.createElement('div');
        card.className = 'chart-card mb-2';
        card.style.animation = `fadeUp 0.5s ease ${{cIdx * 0.05}}s both`;

        const contigLen = contig.length || 1;
        const regions = contig.regions || [];
        const nRegions = regions.length;

        const displayRegions = regions;

        const headerHtml = `
            <div class="chart-header">
                <h2>${{contig.header}}</h2>
                <p>${{nRegions}} target region(s) &bull; Contig length: ${{contigLen.toLocaleString()}} bp</p>
            </div>
        `;

        const plotId = `genome-plot-${{cIdx}}`;
        window.genomePlots.push(plotId);
        
        // Fixed vertical formatting for consistency
        const heightPx = 160;
        card.innerHTML = headerHtml + `<div class="chart-body"><div id="${{plotId}}" style="height:${{heightPx}}px;"></div></div>`;
        container.appendChild(card);

        // Build Plotly shapes and annotations
        const shapes = [];
        const annotations = [];
        const yBase = 0;
        const barH = 0.4;
        const yRange = 2.0;

        // Contig backbone
        shapes.push({{
            type: 'rect',
            x0: 0, x1: contigLen,
            y0: yBase - barH, y1: yBase + barH,
            fillcolor: '#e2e8f0',
            line: {{ color: '#cbd5e1', width: 1 }},
            layer: 'below'
        }});

        // Annotations for contig ends
        annotations.push({{
            x: 0, y: yBase, text: '0', showarrow: false,
            font: {{ size: 9, color: '#94a3b8' }}, yshift: -25, xanchor: 'left'
        }});
        annotations.push({{
            x: contigLen, y: yBase, text: contigLen.toLocaleString() + ' bp', showarrow: false,
            font: {{ size: 9, color: '#94a3b8' }}, yshift: -25, xanchor: 'right'
        }});

        // Use consistent colors for all regions
        const regionColor = 'rgba(37,99,235,0.25)';
        const regionBorder = '{PALETTE["primary"]}';

        const primerTraces = [
            {{ x: [], y: [], hoverinfo: 'skip', type: 'scatter', mode: 'markers', marker: {{ symbol: 'triangle-right', size: 10, color: '{PALETTE["success"]}', line: {{ color: 'white', width: 1 }} }}, showlegend: false }},
            {{ x: [], y: [], hoverinfo: 'skip', type: 'scatter', mode: 'markers', marker: {{ symbol: 'triangle-left', size: 10, color: '{PALETTE["danger"]}', line: {{ color: 'white', width: 1 }} }}, showlegend: false }}
        ];

        displayRegions.forEach((region, rIdx) => {{
            if (region.start === undefined || region.end === undefined) return;
            const rStart = region.start;
            const rEnd = region.end;

            // Region rectangle on the contig bar
            shapes.push({{
                type: 'rect',
                x0: rStart, x1: rEnd,
                y0: yBase - barH, y1: yBase + barH,
                fillcolor: regionColor,
                line: {{ color: regionBorder, width: 1.5 }},
                layer: 'above'
            }});

            // Region label (R1, R2, ...) placed below the bar
            const midpoint = (rStart + rEnd) / 2;
            const regionLen = rEnd - rStart;
            let tooltipText = `Region ${{rIdx + 1}} (${{regionLen.toLocaleString()}} bp)`;
            
            let hasAnnotation = false;
            let geneId = "";
            if (region.annotation && region.annotation !== "no annotation" && region.annotation !== "-" && region.annotation !== "no annotation (file not loaded)") {{
                geneId = region.annotation;
                tooltipText += `<br>Gene ID: ${{geneId}}`;
                hasAnnotation = true;
            }}
            
            annotations.push({{
                x: midpoint,
                y: yBase - barH - 0.2,
                text: `R${{rIdx + 1}}`,
                hovertext: tooltipText,
                showarrow: false,
                font: {{ size: 9, color: regionBorder, family: 'Inter', weight: 'bold' }},
                yshift: 0
            }});

            if (hasAnnotation) {{
                annotations.push({{
                    name: 'gene_annotation',
                    x: midpoint,
                    y: yBase + barH + 0.3,
                    text: geneId,
                    showarrow: false,
                    font: {{ size: 9, color: '{PALETTE["text_muted"]}' }},
                    visible: document.getElementById('toggleAnnotations') ? document.getElementById('toggleAnnotations').checked : true
                }});
            }}

            // Scatter markers for primers (if defined)
            if (region.primer_left_start !== undefined) {{
                primerTraces[0].x.push(rStart + region.primer_left_start);
                primerTraces[0].y.push(yBase);
                
                const rightX = rStart + (region.primer_right_start || regionLen - 20) + (region.primer_right_len || 20);
                primerTraces[1].x.push(rightX);
                primerTraces[1].y.push(yBase);
            }}
        }});

        // Built layout data
        Plotly.newPlot(plotId, [
            {{ x: [null], y: [null], type: 'scatter', mode: 'none', showlegend: false }},
            primerTraces[0],
            primerTraces[1]
        ], {{
            ...layoutBase,
            xaxis: {{
                title: 'Position (bp)',
                range: [-contigLen * 0.02, contigLen * 1.02],
                gridcolor: 'rgba(0,0,0,0.04)',
                zeroline: false,
                tickformat: ',d'
            }},
            yaxis: {{
                range: [-yRange, yRange],
                showticklabels: false,
                showgrid: false,
                zeroline: false
            }},
            shapes: shapes,
            annotations: annotations,
            margin: {{ t: 10, r: 30, l: 30, b: 45 }},
            showlegend: false
        }}, plotConfig);
    }});

    // Add a legend card
    const legendCard = document.createElement('div');
    legendCard.className = 'chart-card';
    legendCard.innerHTML = `
        <div style="padding:1rem 1.5rem; display:flex; flex-wrap:wrap; gap:1.5rem; align-items:center; font-size:0.8rem; color:var(--text-muted);">
            <span><span style="display:inline-block;width:14px;height:14px;background:#e2e8f0;border:1px solid #cbd5e1;border-radius:3px;vertical-align:middle;margin-right:4px;"></span> Contig backbone</span>
            <span><span style="display:inline-block;width:14px;height:14px;background:rgba(37,99,235,0.25);border:1px solid {PALETTE["primary"]};border-radius:3px;vertical-align:middle;margin-right:4px;"></span> Target Region (R1, R2...)</span>
            <span><span style="display:inline-block;width:0;height:0;border-top:5px solid transparent;border-bottom:5px solid transparent;border-left:7px solid {PALETTE["success"]};vertical-align:middle;margin-right:4px;"></span> Forward primer</span>
            <span><span style="display:inline-block;width:0;height:0;border-top:5px solid transparent;border-bottom:5px solid transparent;border-right:7px solid {PALETTE["danger"]};vertical-align:middle;margin-right:4px;"></span> Reverse primer</span>
        </div>
    `;
    container.appendChild(legendCard);
}})();

</script>
</body>
</html>
"""

    report_path = os.path.join(results_dir, "diffprimer_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    logger.info(f"Report generated successfully: [bold white]{report_path}[/bold white]")
