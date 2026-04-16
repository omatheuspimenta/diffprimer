import os
import re
import json
import pandas as pd
from diffprimer.logs import diffprimerLog

logger = diffprimerLog()

def parse_header(header_col):
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

def generate_html_report(results_dir: str):
    csv_file = os.path.join(results_dir, "primer_design_results.csv")
    
    if not os.path.exists(csv_file):
        logger.error(f"Cannot generate report: {csv_file} does not exist.")
        return

    logger.info("Generating HTML report...")
    
    # Load data
    try:
        df = pd.read_csv(csv_file, sep=";")
    except Exception as e:
        logger.error(f"Failed to read CSV file: {e}")
        return

    # Create parsed columns from Header
    parsed_header_df = parse_header(df["Header"])
    
    # Ensure all required columns exist, fill with NA if missing
    expected_cols = [
        "Header", "Forw_Seq", "Rev_Seq", "Amplicon_Seq", "Amplicon_Size",
        "Region_Length", "Annotation", "Sequence_Region", "Specificity_Tag",
        "Most_Similar_Target", "Max_Similarity_pct"
    ]
    for col in expected_cols:
        if col not in df.columns:
            df[col] = "NA"
            
    # Combine parsed columns and expected columns for the table display
    table_df = pd.concat([parsed_header_df, df[expected_cols]], axis=1)
    
    # Pre-compute data for charting
    # Amplicon Sizes (clean NA)
    amp_sizes = df["Amplicon_Size"].replace("NA", pd.NA).dropna().astype(float).tolist()
    
    # Specificity Tags
    spec_tags = df["Specificity_Tag"].value_counts().to_dict()
    spec_labels = list(spec_tags.keys())
    spec_values = list(spec_tags.values())
    
    # TM and GC
    left_tm = df.get("Left_TM", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    right_tm = df.get("Right_TM", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    left_gc = df.get("Left_GC", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    right_gc = df.get("Right_GC", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    
    # Hairpin and Dimerization (Boxplots)
    left_hairpin = df.get("Left_HAIRPIN_TH", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    right_hairpin = df.get("Right_HAIRPIN_TH", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    left_self = df.get("Left_SELF_ANY_TH", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    right_self = df.get("Right_SELF_ANY_TH", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()
    
    # Specificity similarity percentages for Histogram
    max_sim_pct = df.get("Max_Similarity_pct", pd.Series()).replace("NA", pd.NA).dropna().astype(float).tolist()

    # Generate HTML string
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>diffprimer Report</title>
    <!-- Google Fonts -->
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap" rel="stylesheet">
    <!-- DataTables CSS -->
    <link rel="stylesheet" href="https://cdn.datatables.net/2.0.3/css/dataTables.dataTables.css" />
    <style>
        :root {{
            --bg-color: #f4f7f6;
            --surface-color: #ffffff;
            --text-color: #333333;
            --text-muted: #666666;
            --primary-color: #2b6cb0; /* Professional Blue */
            --border-color: #e2e8f0;
            --glass-bg: rgba(255, 255, 255, 0.7);
            --glass-border: rgba(255, 255, 255, 0.5);
            --shadow: 0 4px 6px rgba(0, 0, 0, 0.05), 0 1px 3px rgba(0, 0, 0, 0.1);
        }}

        @media (prefers-color-scheme: dark) {{
            /* Optional basic dark mode support if system set, but design is primarily light */
            :root.dark-theme {{
                --bg-color: #1a202c;
                --surface-color: #2d3748;
                --text-color: #edf2f7;
                --text-muted: #a0aec0;
                --primary-color: #63b3ed;
                --border-color: #4a5568;
                --glass-bg: rgba(45, 55, 72, 0.8);
                --glass-border: rgba(255, 255, 255, 0.1);
                --shadow: 0 4px 6px rgba(0, 0, 0, 0.2);
            }}
        }}

        body {{
            font-family: 'Inter', sans-serif;
            background-color: var(--bg-color);
            color: var(--text-color);
            margin: 0;
            padding: 0;
            transition: background-color 0.3s ease, color 0.3s ease;
        }}

        header {{
            background: linear-gradient(135deg, #2b6cb0 0%, #4299e1 100%);
            color: white;
            padding: 2rem 5%;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            position: relative;
        }}

        .header-content {{
            max-width: 1400px;
            margin: 0 auto;
            position: relative;
            z-index: 2;
        }}

        h1 {{ margin: 0; font-size: 2.5rem; font-weight: 700; letter-spacing:-0.5px; }}
        .header-subtitle {{ margin-top: 0.5rem; font-size: 1.1rem; opacity: 0.9; }}

        .container {{
            max-width: 1400px;
            margin: 2rem auto;
            padding: 0 5%;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}

        .stat-card {{
            background: var(--surface-color);
            background: var(--glass-bg);
            backdrop-filter: blur(10px);
            border: 1px solid var(--glass-border);
            border-radius: 12px;
            padding: 1.5rem;
            box-shadow: var(--shadow);
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }}
        .stat-card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 6px 15px rgba(0, 0, 0, 0.08);
        }}

        .stat-card h3 {{ margin: 0 0 0.5rem 0; font-size: 0.9rem; color: var(--text-muted); text-transform: uppercase; letter-spacing: 1px; }}
        .stat-card .value {{ font-size: 2.2rem; font-weight: 700; color: var(--primary-color); }}

        .charts-container {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 2rem;
            margin-bottom: 3rem;
        }}

        .chart-box {{
            background: var(--surface-color);
            border-radius: 12px;
            padding: 1.5rem;
            box-shadow: var(--shadow);
            border: 1px solid var(--border-color);
        }}
        .chart-box h2 {{ margin-top: 0; font-size: 1.3rem; border-bottom: 1px solid var(--border-color); padding-bottom: 0.5rem; color: var(--text-color); }}

        .table-container {{
            background: var(--surface-color);
            border-radius: 12px;
            padding: 2rem;
            box-shadow: var(--shadow);
            border: 1px solid var(--border-color);
            overflow-x: auto;
        }}
        
        .table-container h2 {{ margin-top: 0; font-size: 1.5rem; border-bottom: 1px solid var(--border-color); padding-bottom: 1rem; color: var(--text-color); margin-bottom: 1.5rem; }}

        /* DataTables Custom Overrides */
        table.dataTable {{
            border-collapse: collapse !important;
            width: 100% !important;
            color: var(--text-color);
        }}
        table.dataTable th, table.dataTable td {{
            font-family: 'Inter', sans-serif;
            font-size: 0.9rem;
            padding: 12px 15px;
            border-bottom: 1px solid var(--border-color);
        }}
        table.dataTable thead th {{
            background-color: rgba(0,0,0,0.02);
            color: var(--text-muted);
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            font-size: 0.8rem;
            border-bottom: 2px solid var(--border-color);
        }}
        table.dataTable tbody tr {{
            background-color: transparent !important;
            transition: background-color 0.15s;
        }}
        table.dataTable tbody tr:hover {{
            background-color: rgba(43, 108, 176, 0.04) !important;
        }}
        .dt-search input {{
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 8px 12px;
            background: var(--surface-color);
            color: var(--text-color);
            font-family: 'Inter', sans-serif;
        }}
        
        .seq-cell {{
            font-family: monospace;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            max-width: 150px;
            display: inline-block;
            vertical-align: middle;
        }}
        
        .copy-btn {{
            background: none;
            border: none;
            cursor: pointer;
            color: var(--text-muted);
            margin-left: 6px;
            padding: 2px;
            opacity: 0.5;
            transition: opacity 0.2s, color 0.2s;
            vertical-align: middle;
            display: inline-flex;
            align-items: center;
        }}
        .copy-btn:hover {{
            opacity: 1;
            color: var(--primary-color);
        }}
        .copy-btn svg {{
            width: 14px;
            height: 14px;
        }}
        
        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 12px;
            font-size: 0.75rem;
            font-weight: 600;
            background: #edf2f7;
            color: #4a5568;
        }}
        .badge.specific {{ background: #def7ec; color: #046c4e; }}
        .badge.nonspecific {{ background: #fde8e8; color: #9b1c1c; }}
        .badge.notchecked {{ background: #edf2f7; color: #4a5568; }}

        /* Utility */
        .flex-row {{ display: flex; align-items: center; justify-content: space-between; }}
    </style>
</head>
<body>

    <header>
        <div class="header-content">
            <h1>diffprimer Report</h1>
            <div class="header-subtitle">Interactive Analysis of Target-Specific PCR Primers</div>
        </div>
    </header>

    <div class="container">
        
        <div class="stats-grid">
            <div class="stat-card">
                <h3>Total Candidate Regions</h3>
                <div class="value">{len(df)}</div>
            </div>
            <div class="stat-card">
                <h3>Designed Pairs</h3>
                <div class="value">{len(amp_sizes)}</div>
            </div>
            <div class="stat-card">
                <h3>Specific Pairs</h3>
                <div class="value">{spec_tags.get("Specific_HighSim", 0) + spec_tags.get("Specific_LowGlobalSim", 0)}</div>
            </div>
        </div>

        <div class="charts-container">
            <div class="chart-box">
                <h2>Amplicon Size Distribution</h2>
                <div id="plot-amplicon"></div>
            </div>
            <div class="chart-box">
                <h2>Specificity Breakdown</h2>
                <div id="plot-specificity"></div>
            </div>
            <div class="chart-box">
                <h2>Max Homology Similarity (%)</h2>
                <div id="plot-homology"></div>
            </div>
            <div class="chart-box">
                <h2>Dimerization Risks</h2>
                <div id="plot-dimer"></div>
            </div>
            <div class="chart-box" style="grid-column: 1 / -1; max-width: 1000px; justify-self: center; width: 100%;">
                <h2>Primer Parameter Space (Tm vs GC)</h2>
                <div id="plot-tmgc"></div>
            </div>
        </div>

        <div class="table-container">
            <h2>Detailed Primer Data</h2>
            <table id="primerTable" class="display" style="width:100%">
                <thead>
                    <tr>
                        <th>Sequence ID</th>
                        <th>Region ID</th>
                        <th>Start</th>
                        <th>End</th>
                        <th>Header</th>
                        <th>Forw Seq</th>
                        <th>Rev Seq</th>
                        <th>Amplicon Seq</th>
                        <th>Amplicon Size</th>
                        <th>Region Length</th>
                        <th>Annotation</th>
                        <th>Sequence Region</th>
                        <th>Specificity Tag</th>
                        <th>Most Similar Target</th>
                        <th>Max Sim%</th>
                    </tr>
                </thead>
                <tbody>
"""

    def copy_cell(full_value, display_value, truncate=False):
        if truncate:
            display_html = f'<span class="seq-cell" title="{full_value}">{display_value}</span>'
        else:
            display_html = f'<span>{display_value}</span>'
            
        full_value_escaped = str(full_value).replace("'", "\\'")
        copy_icon = '''<svg fill="currentColor" viewBox="0 0 20 20"><path d="M8 3a1 1 0 011-1h2a1 1 0 110 2H9a1 1 0 01-1-1z"></path><path d="M6 3a2 2 0 00-2 2v11a2 2 0 002 2h8a2 2 0 002-2V5a2 2 0 00-2-2 3 3 0 01-3 3H9a3 3 0 01-3-3z"></path></svg>'''
        
        return f'{display_html} <button class="copy-btn" onclick="copyText(this, \'{full_value_escaped}\')" title="Copy"> {copy_icon} </button>'

    # Populate rows
    for i, row in table_df.iterrows():
        tag = row.get("Specificity_Tag", "Not_Checked")
        badge_class = "notchecked"
        if "Specific" in str(tag):
            badge_class = "specific"
        elif "NonSpecific" in str(tag):
            badge_class = "nonspecific"
            
        # Format the specific badge with a copy button
        tag_html = f'<span class="badge {badge_class}">{tag}</span> <button class="copy-btn" onclick="copyText(this, \'{tag}\')" title="Copy"><svg fill="currentColor" viewBox="0 0 20 20"><path d="M8 3a1 1 0 011-1h2a1 1 0 110 2H9a1 1 0 01-1-1z"></path><path d="M6 3a2 2 0 00-2 2v11a2 2 0 002 2h8a2 2 0 002-2V5a2 2 0 00-2-2 3 3 0 01-3 3H9a3 3 0 01-3-3z"></path></svg></button>'
            
        html_content += f"""
                    <tr>
                        <td>{copy_cell(row.get("Sequence_ID", ""), row.get("Sequence_ID", ""))}</td>
                        <td>{copy_cell(row.get("Region_ID", ""), row.get("Region_ID", ""))}</td>
                        <td>{copy_cell(row.get("Start", ""), row.get("Start", ""))}</td>
                        <td>{copy_cell(row.get("End", ""), row.get("End", ""))}</td>
                        <td>{copy_cell(row.get("Header", ""), row.get("Header", ""), truncate=True)}</td>
                        <td>{copy_cell(row.get("Forw_Seq", ""), row.get("Forw_Seq", ""), truncate=True)}</td>
                        <td>{copy_cell(row.get("Rev_Seq", ""), row.get("Rev_Seq", ""), truncate=True)}</td>
                        <td>{copy_cell(row.get("Amplicon_Seq", ""), row.get("Amplicon_Seq", ""), truncate=True)}</td>
                        <td>{copy_cell(row.get("Amplicon_Size", ""), row.get("Amplicon_Size", ""))}</td>
                        <td>{copy_cell(row.get("Region_Length", ""), row.get("Region_Length", ""))}</td>
                        <td>{copy_cell(row.get("Annotation", ""), row.get("Annotation", ""), truncate=True)}</td>
                        <td>{copy_cell(row.get("Sequence_Region", ""), row.get("Sequence_Region", ""), truncate=True)}</td>
                        <td>{tag_html}</td>
                        <td>{copy_cell(row.get("Most_Similar_Target", ""), row.get("Most_Similar_Target", ""), truncate=True)}</td>
                        <td>{copy_cell(row.get("Max_Similarity_pct", ""), row.get("Max_Similarity_pct", ""))}</td>
                    </tr>
"""

    html_content += f"""
                </tbody>
            </table>
        </div>

    </div>

    <!-- Scripts -->
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    <script src="https://cdn.datatables.net/2.0.3/js/dataTables.js"></script>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    
    <script>
        // Copy to clipboard functionality
        function copyText(button, text) {{
            navigator.clipboard.writeText(text).then(() => {{
                const originalSvg = button.innerHTML;
                // Checkmark svg
                button.innerHTML = '<svg fill="currentColor" viewBox="0 0 20 20" style="color: #046c4e;"><path fill-rule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clip-rule="evenodd"></path></svg>';
                setTimeout(() => {{
                    button.innerHTML = originalSvg;
                }}, 2000);
            }}).catch(err => {{
                console.error('Failed to copy: ', err);
            }});
        }}

        $(document).ready( function () {{
            $('#primerTable').DataTable({{
                pageLength: 25,
                scrollX: true,
                language: {{
                    search: "Filter records:"
                }}
            }});
        }} );

        // Shared Plotly Config for Aesthetics
        const layoutConfig = {{
            font: {{ family: 'Inter, sans-serif' }},
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            margin: {{ t: 20, r: 20, l: 50, b: 50 }}
        }};

        // Plot 1: Amplicon Histogram
        var trace_amp = {{
            x: {json.dumps(amp_sizes)},
            type: 'histogram',
            marker: {{ color: '#3182ce', line: {{ color: '#2b6cb0', width: 1}} }},
            opacity: 0.8
        }};
        Plotly.newPlot('plot-amplicon', [trace_amp], {{
            ...layoutConfig,
            xaxis: {{ title: 'Amplicon Size (bp)' }},
            yaxis: {{ title: 'Count' }}
        }});

        // Plot 2: Specificity Donut
        var trace_spec = {{
            labels: {json.dumps(spec_labels)},
            values: {json.dumps(spec_values)},
            type: 'pie',
            hole: 0.5,
            marker: {{ colors: ['#def7ec', '#fde8e8', '#edf2f7', '#ebf8ff', '#f0fff4', '#fff5f5'] }},
            textinfo: 'label+percent'
        }};
        Plotly.newPlot('plot-specificity', [trace_spec], {{
            ...layoutConfig,
            showlegend: false,
            margin: {{ t: 10, b: 10, l: 10, r: 10 }}
        }});

        // Plot 3: Tm vs GC Scatter
        var trace_tm_gc_left = {{
            x: {json.dumps(left_gc)},
            y: {json.dumps(left_tm)},
            mode: 'markers',
            type: 'scatter',
            name: 'Forward Primer',
            marker: {{ size: 8, color: '#3182ce', opacity: 0.7 }}
        }};
        var trace_tm_gc_right = {{
            x: {json.dumps(right_gc)},
            y: {json.dumps(right_tm)},
            mode: 'markers',
            type: 'scatter',
            name: 'Reverse Primer',
            marker: {{ size: 8, color: '#dd6b20', opacity: 0.7 }}
        }};
        Plotly.newPlot('plot-tmgc', [trace_tm_gc_left, trace_tm_gc_right], {{
            ...layoutConfig,
            xaxis: {{ title: 'GC Content (%)' }},
            yaxis: {{ title: 'Melting Temperature (Tm °C)' }},
            legend: {{ orientation: 'h', y: -0.2 }}
        }});

        // Plot 4: Max Similarity Homology Histogram
        var trace_homology = {{
            x: {json.dumps(max_sim_pct)},
            type: 'histogram',
            marker: {{ color: '#805ad5', line: {{ color: '#553c9a', width: 1 }} }},
            opacity: 0.8
        }};
        Plotly.newPlot('plot-homology', [trace_homology], {{
            ...layoutConfig,
            xaxis: {{ title: 'Similarity to Off-Target (%)' }},
            yaxis: {{ title: 'Count' }}
        }});

        // Plot 5: Dimerization & Hairpin Boxplot
        var yData_dimer = {json.dumps(left_hairpin)}
            .map(x => "F-Hairpin")
            .concat({json.dumps(right_hairpin)}.map(x => "R-Hairpin"))
            .concat({json.dumps(left_self)}.map(x => "F-Self"))
            .concat({json.dumps(right_self)}.map(x => "R-Self"));
            
        var xData_dimer = {json.dumps(left_hairpin)}
            .concat({json.dumps(right_hairpin)})
            .concat({json.dumps(left_self)})
            .concat({json.dumps(right_self)});

        var trace_dimer = {{
            x: xData_dimer,
            y: yData_dimer,
            type: 'box',
            orientation: 'h',
            marker: {{ color: '#4fd1c5' }},
            line: {{ width: 1 }}
        }};
        Plotly.newPlot('plot-dimer', [trace_dimer], {{
            ...layoutConfig,
            xaxis: {{ title: 'Propensity Score' }},
            margin: {{ t: 20, r: 20, b: 50, l: 80 }}
        }});

    </script>
</body>
</html>
"""

    report_path = os.path.join(results_dir, "diffprimer_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html_content)
        
    logger.info(f"Report generated successfully: [bold white]{report_path}[/bold white]")
