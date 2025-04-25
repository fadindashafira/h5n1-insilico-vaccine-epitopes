import streamlit as st
import pandas as pd
import numpy as np
import os
import subprocess
import time
import glob
import re
import base64
from datetime import datetime
import plotly.express as px
import plotly.graph_objects as go
from PIL import Image
import json
import threading
import queue

# Set page configuration with Kalbe theme
st.set_page_config(
    page_title="Kalbe Bioinformatics - H5N1 Vaccine Design",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Add custom CSS for Kalbe theme (green)
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2e7d32;  /* Kalbe green */
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.5rem;
        font-weight: 600;
        color: #43a047;  /* Lighter green */
        margin-bottom: 0.5rem;
    }
    .info-box {
        background-color: #e8f5e9;  /* Very light green */
        border-left: 5px solid #43a047;
        padding: 1rem;
        margin-bottom: 1rem;
    }
    .success-box {
        background-color: #e8f5e9;
        border-left: 5px solid #2e7d32;
        padding: 1rem;
        margin-bottom: 1rem;
    }
    .warning-box {
        background-color: #fff8e1;
        border-left: 5px solid #ffb300;
        padding: 1rem;
        margin-bottom: 1rem;
    }
    .error-box {
        background-color: #ffebee;
        border-left: 5px solid #c62828;
        padding: 1rem;
        margin-bottom: 1rem;
    }
    .terminal {
        background-color: #212121;
        color: #f5f5f5;
        padding: 1rem;
        border-radius: 0.5rem;
        font-family: monospace;
        overflow-x: auto;
        margin-bottom: 1rem;
    }
    .stProgress > div > div > div > div {
        background-color: #43a047;  /* Progress bar color - green */
    }
    .stButton>button {
        background-color: #43a047;
        color: white;
    }
    .stButton>button:hover {
        background-color: #2e7d32;
        color: white;
    }
    .sidebar .sidebar-content {
        background-color: #e8f5e9;
    }
    /* Custom metrics styling */
    [data-testid="stMetricValue"] {
        color: #2e7d32;
        font-weight: bold;
    }
    /* Custom sidebar styling */
    [data-testid="stSidebar"] {
        background-color: #f1f8e9;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state for storing results without needing to refresh
if 'results' not in st.session_state:
    st.session_state.results = {
        'bcell_epitopes': [],
        'tcell_i_epitopes': [],
        'tcell_ii_epitopes': [],
        'combined_epitopes': [],
        'vaccine_constructs': []
    }

if 'log_output' not in st.session_state:
    st.session_state.log_output = ""

if 'pipeline_status' not in st.session_state:
    st.session_state.pipeline_status = "Not Started"

if 'progress' not in st.session_state:
    st.session_state.progress = 0.0

if 'current_task' not in st.session_state:
    st.session_state.current_task = "No task running"

# Helper functions
def load_kalbe_logo():
    """Load and encode the Kalbe Bioinformatics logo"""
    logo_path = "./img/kalbe_bioinformatics_logo.jpeg"
    if os.path.exists(logo_path):
        return Image.open(logo_path)
    else:
        # Create a placeholder image if the logo file doesn't exist
        return None

def run_pipeline_async(cmd, output_queue):
    """Run pipeline asynchronously and put results in queue"""
    try:
        st.session_state.pipeline_status = "Running"
        st.session_state.log_output = ""
        st.session_state.progress = 0.0
        st.session_state.current_task = "Starting pipeline..."
        
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            text=True
        )
        
        # Monitor the process output
        while True:
            line = process.stdout.readline()
            if not line and process.poll() is not None:
                break
            if line:
                st.session_state.log_output += line
                
                # Update progress based on the output
                if "Process complete" in line:
                    progress_match = re.search(r'(\d+)%', line)
                    if progress_match:
                        progress = float(progress_match.group(1)) / 100.0
                        st.session_state.progress = progress
                
                # Update current task
                task_match = re.search(r'Running: ([^\\n]+)', line)
                if task_match:
                    st.session_state.current_task = task_match.group(1)
        
        # Check if process completed successfully
        if process.returncode == 0:
            st.session_state.pipeline_status = "Completed"
            st.session_state.progress = 1.0
            
            # Load the results directly after completion
            experiment_id = cmd.split("--experiment_id")[1].split()[0].strip()
            experiment_output = f"results/results_{experiment_id}"
            results = get_latest_results(experiment_output)
            output_queue.put({"status": "success", "results": results})
        else:
            error_output = process.stderr.read()
            st.session_state.pipeline_status = "Failed"
            output_queue.put({"status": "error", "message": error_output})
    
    except Exception as e:
        st.session_state.pipeline_status = "Error"
        output_queue.put({"status": "error", "message": str(e)})

def parse_bcpreds_output(bcpreds_output):
    """Parse BCPREDS output and convert to DataFrame"""
    epitopes = []
    
    lines = bcpreds_output.strip().split('\n')
    for line in lines:
        if line.strip() and not line.startswith("Position") and not line.startswith("-"):
            parts = line.strip().split()
            if len(parts) >= 3:
                position = int(parts[0])
                epitope = parts[1]
                score = float(parts[2])
                epitopes.append({
                    "Position": position,
                    "Epitope": epitope,
                    "Score": score
                })
    
    return pd.DataFrame(epitopes)

def get_latest_results(experiment_output):
    """Get the latest results files from the experiment output directory"""
    if not os.path.exists(experiment_output):
        # For demo purposes, load sample results if directory doesn't exist
        return load_sample_results()
    
    results = {
        'bcell_epitopes': [],
        'tcell_i_epitopes': [],
        'tcell_ii_epitopes': [],
        'combined_epitopes': [],
        'vaccine_constructs': []
    }
    
    # Find B-cell epitope files
    bcell_files = glob.glob(f"{experiment_output}/*_bcell_epitopes.csv")
    for file in bcell_files:
        try:
            protein = os.path.basename(file).split('_')[0]
            df = pd.read_csv(file)
            results['bcell_epitopes'].append({
                'protein': protein,
                'file': file,
                'data': df
            })
        except:
            pass
    
    # Find T-cell Class I epitope files
    tcell_i_files = glob.glob(f"{experiment_output}/*_tcell_i_epitopes.csv")
    for file in tcell_i_files:
        try:
            protein = os.path.basename(file).split('_')[0]
            df = pd.read_csv(file)
            results['tcell_i_epitopes'].append({
                'protein': protein,
                'file': file,
                'data': df
            })
        except:
            pass
    
    # Find T-cell Class II epitope files
    tcell_ii_files = glob.glob(f"{experiment_output}/*_tcell_ii_epitopes.csv")
    for file in tcell_ii_files:
        try:
            protein = os.path.basename(file).split('_')[0]
            df = pd.read_csv(file)
            results['tcell_ii_epitopes'].append({
                'protein': protein,
                'file': file,
                'data': df
            })
        except:
            pass
    
    # Find combined epitope files
    combined_files = glob.glob(f"{experiment_output}/*_combined_epitopes.csv")
    for file in combined_files:
        try:
            protein = os.path.basename(file).split('_')[0]
            df = pd.read_csv(file)
            results['combined_epitopes'].append({
                'protein': protein,
                'file': file,
                'data': df
            })
        except:
            pass
    
    # Find vaccine construct files
    vaccine_files = glob.glob(f"{experiment_output}/*_vaccine_construct.fasta")
    for file in vaccine_files:
        try:
            protein = os.path.basename(file).split('_')[0]
            with open(file, 'r') as f:
                content = f.read()
            results['vaccine_constructs'].append({
                'protein': protein,
                'file': file,
                'content': content
            })
        except:
            pass
    
    return results

def load_sample_results():
    """Load sample results for demo purposes"""
    # Sample B-cell epitopes for hemagglutinin based on paper
    ha_bcell = pd.DataFrame([
        {"Position": 98, "Epitope": "KANPNNDLC", "Score": 0.998},
        {"Position": 137, "Epitope": "SWSDHEASS", "Score": 0.997},
        {"Position": 181, "Epitope": "NNTNQEDLL", "Score": 0.996},
        {"Position": 338, "Epitope": "QRESRRKKR", "Score": 0.985},
        {"Position": 216, "Epitope": "IGTSTLNQR", "Score": 0.984},
        {"Position": 288, "Epitope": "GNCNTKCQT", "Score": 0.979},
        {"Position": 510, "Epitope": "EEARLKREE", "Score": 0.973},
        {"Position": 10, "Epitope": "MVSLVKSDQ", "Score": 0.964},
        {"Position": 557, "Epitope": "CSNGSLQCR", "Score": 0.935}
    ])
    
    # Sample T-cell Class I epitopes for hemagglutinin
    ha_tcell_i = pd.DataFrame([
        {"Position": 1, "Epitope": "MEKIVLLLA", "Allele": "HLA-B*35:01", "Score": 0.85, "Affinity": "Strong Binder"},
        {"Position": 2, "Epitope": "EKIVLLLAM", "Allele": "HLA-B*35:01", "Score": 0.75, "Affinity": "Weak Binder"},
        {"Position": 151, "Epitope": "CPYLGSPSF", "Allele": "HLA-B*07:02", "Score": 0.92, "Affinity": "Strong Binder"},
        {"Position": 293, "Epitope": "KCQTPMGAI", "Allele": "HLA-B*07:02", "Score": 0.88, "Affinity": "Strong Binder"},
        {"Position": 389, "Epitope": "KAVDGVTNK", "Allele": "HLA-A*11:01", "Score": 0.77, "Affinity": "Weak Binder"}
    ])
    
    # Sample T-cell Class II epitopes for hemagglutinin
    ha_tcell_ii = pd.DataFrame([
        {"Position": 10, "Epitope": "MVSLVKSDQ", "Allele": "HLA-DRB1*03:01", "Score": 0.65, "Affinity": "Weak Binder"},
        {"Position": 216, "Epitope": "IGTSTLNQR", "Allele": "HLA-DRB1*03:01", "Score": 0.81, "Affinity": "Strong Binder"}
    ])
    
    # Sample B-cell epitopes for neuraminidase
    na_bcell = pd.DataFrame([
        {"Position": 342, "Epitope": "TKSTNSRSG", "Score": 0.996},
        {"Position": 174, "Epitope": "GISGPDNEA", "Score": 0.995},
        {"Position": 142, "Epitope": "PVGEAPSPY", "Score": 0.995},
        {"Position": 303, "Epitope": "GDNPRPNDG", "Score": 0.989},
        {"Position": 84, "Epitope": "NNIRIGSKG", "Score": 0.977},
        {"Position": 188, "Epitope": "YNGIITDTI", "Score": 0.975},
        {"Position": 257, "Epitope": "EESCSCYDA", "Score": 0.972}
    ])
    
    # Sample T-cell Class I epitopes for neuraminidase
    na_tcell_i = pd.DataFrame([
        {"Position": 2, "Epitope": "NPNQKIITI", "Allele": "HLA-B*07:02", "Score": 0.89, "Affinity": "Strong Binder"},
        {"Position": 261, "Epitope": "CYPDAGEIT", "Allele": "HLA-A*24", "Score": 0.79, "Affinity": "Strong Binder"},
        {"Position": 398, "Epitope": "IRPCFWVEL", "Allele": "HLA-B*27:05", "Score": 0.86, "Affinity": "Strong Binder"},
        {"Position": 399, "Epitope": "RPCFWVELI", "Allele": "HLA-B*07:02", "Score": 0.72, "Affinity": "Weak Binder"}
    ])
    
    # Sample T-cell Class II epitopes for neuraminidase
    na_tcell_ii = pd.DataFrame([
        {"Position": 188, "Epitope": "YNGIITDTI", "Allele": "HLA-DRB1*04:05", "Score": 0.68, "Affinity": "Weak Binder"}
    ])
    
    # Sample combined epitopes for hemagglutinin
    ha_combined = pd.DataFrame([
        {"Position": 216, "Epitope": "IGTSTLNQR", "Type": "B-cell & T-cell II", "Score": 0.984, "HLA": "HLA-DRB1*03:01", "ΔG binding": -56.958},
        {"Position": 10, "Epitope": "MVSLVKSDQ", "Type": "B-cell & T-cell II", "Score": 0.964, "HLA": "HLA-DRB1*03:01", "ΔG binding": -11.776},
        {"Position": 1, "Epitope": "MEKIVLLLA", "Type": "T-cell I", "Score": 0.85, "HLA": "HLA-B*35:01", "ΔG binding": -40.015},
        {"Position": 151, "Epitope": "CPYLGSPSF", "Type": "T-cell I", "Score": 0.92, "HLA": "HLA-B*07:02", "ΔG binding": -28.047},
        {"Position": 293, "Epitope": "KCQTPMGAI", "Type": "T-cell I", "Score": 0.88, "HLA": "HLA-B*07:02", "ΔG binding": -33.611}
    ])
    
    # Sample combined epitopes for neuraminidase
    na_combined = pd.DataFrame([
        {"Position": 188, "Epitope": "YNGIITDTI", "Type": "B-cell & T-cell II", "Score": 0.975, "HLA": "HLA-DRB1*04:05", "ΔG binding": -17.805},
        {"Position": 2, "Epitope": "NPNQKIITI", "Type": "T-cell I", "Score": 0.89, "HLA": "HLA-B*07:02", "ΔG binding": -23.878},
        {"Position": 261, "Epitope": "CYPDAGEIT", "Type": "T-cell I", "Score": 0.79, "HLA": "HLA-A*24", "ΔG binding": -33.020}
    ])
    
    # Sample vaccine construct for hemagglutinin
    ha_vaccine = """>Hemagglutinin_Epitope_Vaccine_H5N1_Indonesian_Strain
IGTSTLNQRGPGPGMVSLVKSDQGPGPGMEKIVLLLACPYLGSPSFGPGPGKCQTPMGAI"""
    
    # Sample vaccine construct for neuraminidase
    na_vaccine = """>Neuraminidase_Epitope_Vaccine_H5N1_Indonesian_Strain
YNGIITDTIGPGPGNPNQKIITIGPGPGCYPDAGEIT"""
    
    # Sample vaccine construct for combined
    combined_vaccine = """>Combined_Epitope_Vaccine_H5N1_Indonesian_Strain
IGTSTLNQRGPGPGMVSLVKSDQGPGPGMEKIVLLLACPYLGSPSFGPGPGKCQTPMGAIGPGPGYNGIITDTIGPGPGNPNQKIITIGPGPGCYPDAGEIT"""
    
    return {
        'bcell_epitopes': [
            {'protein': 'hemagglutinin', 'file': 'ha_bcell.csv', 'data': ha_bcell},
            {'protein': 'neuraminidase', 'file': 'na_bcell.csv', 'data': na_bcell}
        ],
        'tcell_i_epitopes': [
            {'protein': 'hemagglutinin', 'file': 'ha_tcell_i.csv', 'data': ha_tcell_i},
            {'protein': 'neuraminidase', 'file': 'na_tcell_i.csv', 'data': na_tcell_i}
        ],
        'tcell_ii_epitopes': [
            {'protein': 'hemagglutinin', 'file': 'ha_tcell_ii.csv', 'data': ha_tcell_ii},
            {'protein': 'neuraminidase', 'file': 'na_tcell_ii.csv', 'data': na_tcell_ii}
        ],
        'combined_epitopes': [
            {'protein': 'hemagglutinin', 'file': 'ha_combined.csv', 'data': ha_combined},
            {'protein': 'neuraminidase', 'file': 'na_combined.csv', 'data': na_combined}
        ],
        'vaccine_constructs': [
            {'protein': 'hemagglutinin', 'file': 'ha_vaccine.fasta', 'content': ha_vaccine},
            {'protein': 'neuraminidase', 'file': 'na_vaccine.fasta', 'content': na_vaccine},
            {'protein': 'combined', 'file': 'combined_vaccine.fasta', 'content': combined_vaccine}
        ]
    }

def main():
    # Load Kalbe logo
    logo = load_kalbe_logo()
    
    # Show logo in the sidebar
    if logo:
        st.sidebar.image(logo, width=250)
    else:
        st.sidebar.title("Kalbe Bioinformatics")
    
    # Main header
    st.markdown('<h1 class="main-header">H5N1 Epitope-Based Vaccine Design</h1>', unsafe_allow_html=True)
    
    # Sidebar for navigation
    st.sidebar.title("Navigation")
    pages = [
        "Home",
        "Configure Pipeline",
        "Run Pipeline",
        "View Results",
        "About"
    ]
    selected_page = st.sidebar.radio("Go to", pages)
    
    # Initialize session state variables if they don't exist
    if 'config' not in st.session_state:
        st.session_state.config = {
            'ha_accession': 'BAL61222.1',  # H5N1 Indonesian strain HA
            'na_accession': 'BAL61230.1',  # H5N1 Indonesian strain NA
            'outdir': 'results',
            'experiment_id': f'exp_{datetime.now().strftime("%Y%m%d_%H%M%S")}',
            'run_md': False
        }
    
    # Check for results from async operations
    if 'result_queue' in st.session_state and not st.session_state.result_queue.empty():
        try:
            result = st.session_state.result_queue.get_nowait()
            if result['status'] == 'success':
                st.session_state.results = result['results']
                st.success("Pipeline completed successfully! Results are ready to view.")
            else:
                st.error(f"Pipeline error: {result['message']}")
        except queue.Empty:
            pass
    
    # Render the selected page
    if selected_page == "Home":
        show_home_page()
    elif selected_page == "Configure Pipeline":
        show_config_page()
    elif selected_page == "Run Pipeline":
        show_run_page()
    elif selected_page == "View Results":
        show_results_page()
    elif selected_page == "About":
        show_about_page()

def show_home_page():
    st.markdown('<h2 class="sub-header">Welcome to the H5N1 Vaccine Design Tool</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
    This application provides a user-friendly interface for designing epitope-based vaccines 
    against H5N1 avian influenza virus using immunoinformatics approaches. The tool 
    automates the prediction of B-cell and T-cell epitopes, analyzes their conservation, 
    and designs optimal vaccine constructs based on methodology from Tambunan et al. (2016).
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    ### Workflow Overview
    
    1. **Input**: Hemagglutinin (HA) and Neuraminidase (NA) protein sequences from H5N1 virus
    2. **Epitope Prediction**: 
       - B-cell epitope prediction using BCPREDS
       - T-cell Class I epitope prediction (ProPred I, NetMHCpan)
       - T-cell Class II epitope prediction (ProPred, NetMHCIIpan)
    3. **Epitope Selection**: Filter epitopes based on binding affinity (ΔG binding)
    4. **Vaccine Design**: Combine selected epitopes with GPGPG linkers
    5. **Molecular Dynamics**: Optional simulation of epitope-HLA interactions
    """)
    
    # Display some key metrics
    st.markdown("### Key Features")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Top Epitopes", "6+", "From the paper")
    
    with col2:
        st.metric("HLA Coverage", "10+ alleles", "Class I & II")
    
    with col3:
        st.metric("Based on", "Indonesian strain", "H5N1 virus")
    
    # Paper reference
    st.markdown("### Reference Publication")
    st.markdown("""
    This tool implements the methodology described in:
    
    **Tambunan et al. (2016). Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions.** *Bioinformatics and Biology Insights*, 10, 27-35.
    
    Key findings from the paper:
    - Identified 4 hemagglutinin epitopes and 2 neuraminidase epitopes with high binding affinity
    - Best hemagglutinin epitopes: MEKIVLLLA, CPYLGSPSF, KCQTPMGAI, IGTSTLNQR
    - Best neuraminidase epitopes: NPNQKIITI, CYPDAGEIT
    """)
    
    # Show a flowchart of the pipeline
    st.markdown("### Pipeline Workflow")
    
    # Create a simple flowchart using Plotly
    fig = go.Figure()
    
    # Add nodes
    nodes = [
        {"name": "Input Sequences", "x": 0, "y": 0},
        {"name": "B-cell Epitope\nPrediction", "x": -1, "y": -1},
        {"name": "T-cell Class I\nEpitope Prediction", "x": 0, "y": -1},
        {"name": "T-cell Class II\nEpitope Prediction", "x": 1, "y": -1},
        {"name": "Epitope Selection\n& Filtering", "x": 0, "y": -2},
        {"name": "Vaccine Design", "x": 0, "y": -3},
        {"name": "Molecular Dynamics\n(Optional)", "x": 0, "y": -4}
    ]
    
    # Add arrows
    arrows = [
        {"from": 0, "to": 1},
        {"from": 0, "to": 2},
        {"from": 0, "to": 3},
        {"from": 1, "to": 4},
        {"from": 2, "to": 4},
        {"from": 3, "to": 4},
        {"from": 4, "to": 5},
        {"from": 5, "to": 6}
    ]
    
    # Plot nodes
    for i, node in enumerate(nodes):
        fig.add_trace(go.Scatter(
            x=[node["x"]],
            y=[node["y"]],
            mode="markers+text",
            marker=dict(size=30, color="#43a047"),
            text=[node["name"]],
            textposition="middle center",
            hoverinfo="text",
            name=node["name"]
        ))
    
    # Plot arrows
    for arrow in arrows:
        fig.add_trace(go.Scatter(
            x=[nodes[arrow["from"]]["x"], nodes[arrow["to"]]["x"]],
            y=[nodes[arrow["from"]]["y"], nodes[arrow["to"]]["y"]],
            mode="lines",
            line=dict(width=2, color="#2e7d32"),
            hoverinfo="none",
            showlegend=False
        ))
    
    fig.update_layout(
        showlegend=False,
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[-1.5, 1.5]
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[-4.5, 0.5]
        ),
        margin=dict(l=20, r=20, t=20, b=20),
        height=500,
        plot_bgcolor="rgba(0,0,0,0)"
    )
    
    st.plotly_chart(fig, use_container_width=True)

def show_config_page():
    st.markdown('<h2 class="sub-header">Configure Pipeline</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
    Configure the parameters for the H5N1 vaccine design pipeline. These settings will be used 
    when executing the pipeline.
    </div>
    """, unsafe_allow_html=True)
    
    # Create tabs for different configuration sections
    tab1, tab2, tab3 = st.tabs(["Input Sequences", "Pipeline Settings", "Advanced Options"])
    
    with tab1:
        st.markdown("### Input Sequences")
        st.markdown("Specify accession numbers for Hemagglutinin and Neuraminidase proteins.")
        
        ha_accession = st.text_input(
            "Hemagglutinin Accession",
            st.session_state.config['ha_accession'],
            help="NCBI accession number for Hemagglutinin protein (e.g., BAL61222.1 for Indonesian H5N1 strain)"
        )
        
        na_accession = st.text_input(
            "Neuraminidase Accession",
            st.session_state.config['na_accession'],
            help="NCBI accession number for Neuraminidase protein (e.g., BAL61230.1 for Indonesian H5N1 strain)"
        )
        
        # Save input sequence settings
        if st.button("Save Sequence Settings"):
            st.session_state.config['ha_accession'] = ha_accession
            st.session_state.config['na_accession'] = na_accession
            st.success("Sequence settings saved successfully!")
    
    with tab2:
        st.markdown("### Pipeline Settings")
        
        # Output directories
        outdir = st.text_input(
            "Output Directory",
            st.session_state.config['outdir'],
            help="Base directory for pipeline output"
        )
        
        experiment_id = st.text_input(
            "Experiment ID",
            st.session_state.config['experiment_id'],
            help="Unique identifier for this experiment run"
        )
        
        # Run molecular dynamics option
        run_md = st.checkbox(
            "Run Molecular Dynamics",
            st.session_state.config.get('run_md', False),
            help="Enable molecular dynamics simulation (requires more computational resources)"
        )
        
        # Save pipeline settings
        if st.button("Save Pipeline Settings"):
            st.session_state.config['outdir'] = outdir
            st.session_state.config['experiment_id'] = experiment_id
            st.session_state.config['run_md'] = run_md
            st.success("Pipeline settings saved successfully!")
    
    with tab3:
        st.markdown("### Advanced Options")
        
        # B-cell epitope prediction
        st.subheader("B-cell Epitope Prediction")
        
        bcell_method_options = ["BCPREDS", "Bepipred-2.0", "Chou-Fasman", "Emini"]
        bcell_method = st.selectbox(
            "Prediction Method",
            bcell_method_options,
            index=bcell_method_options.index("BCPREDS"),
            help="Algorithm for B-cell epitope prediction (BCPREDS used in the paper)"
        )
        
        bcell_threshold = st.slider(
            "Score Threshold",
            0.1, 1.0, 0.8, 0.05,
            help="Minimum score threshold for B-cell epitope prediction (0.8 used in the paper)"
        )
        
        bcell_length = st.slider(
            "Epitope Length",
            8, 20, 9,
            help="Length of B-cell epitopes to predict (9 used in the paper)"
        )
        
        # T-cell epitope prediction
        st.subheader("T-cell Epitope Prediction")
        
        # Class I HLA
        mhci_alleles = st.text_area(
            "Class I HLA Alleles",
            "HLA-A*11:01, HLA-B*07:02, HLA-B*35:01, HLA-B*14:02, HLA-B*39:01",
            help="Comma-separated list of Class I HLA alleles to use for prediction (from paper)"
        )
        
        mhci_methods = st.multiselect(
            "Class I Prediction Methods",
            ["ProPred I", "NetMHCpan"],
            default=["ProPred I", "NetMHCpan"],
            help="Methods used for MHC Class I epitope prediction (both used in paper)"
        )
        
        # Class II HLA
        mhcii_alleles = st.text_area(
            "Class II HLA Alleles",
            "HLA-DRB1*03:01, HLA-DRB1*04:05, HLA-DRB1*12:02, HLA-DRB1*01:01",
            help="Comma-separated list of Class II HLA alleles to use for prediction (from paper)"
        )
        
        mhcii_methods = st.multiselect(
            "Class II Prediction Methods",
            ["ProPred", "NetMHCIIpan"],
            default=["ProPred", "NetMHCIIpan"],
            help="Methods used for MHC Class II epitope prediction (both used in paper)"
        )
        
        # Epitope selection parameters
        st.subheader("Epitope Selection Parameters")
        
        binding_threshold = st.slider(
            "Binding Energy Threshold (kcal/mol)",
            -60.0, 0.0, -20.0, 5.0,
            help="Maximum binding energy (ΔG) threshold for epitope selection (more negative is stronger binding)"
        )
        
        # Vaccine design parameters
        st.subheader("Vaccine Design Parameters")
        
        linker = st.text_input(
            "Epitope Linker",
            "GPGPG",
            help="Amino acid sequence to use as linker between epitopes"
        )
        
        max_epitopes = st.slider(
            "Maximum Epitopes",
            3, 15, 6,
            help="Maximum number of epitopes to include in the vaccine construct"
        )
        
        # These advanced options would be passed to the pipeline
        st.info("These advanced options configure the epitope prediction and selection criteria. Default values are based on the Tambunan et al. (2016) paper.")

def show_run_page():
    st.markdown('<h2 class="sub-header">Run Pipeline</h2>', unsafe_allow_html=True)
    
    # Check if a pipeline is already running
    if st.session_state.pipeline_status == "Running":
        st.markdown("""
        <div class="warning-box">
        A pipeline execution is already in progress. Please wait for it to complete.
        </div>
        """, unsafe_allow_html=True)
        
        # Show progress
        st.subheader("Pipeline Progress")
        
        # Progress bar
        st.progress(st.session_state.progress)
        
        # Current task
        st.info(f"Current task: {st.session_state.current_task}")
        
        # Log output
        with st.expander("Log Output", expanded=True):
            st.markdown(f'<div class="terminal">{st.session_state.log_output}</div>', unsafe_allow_html=True)
        
        # Option to stop the pipeline
        if st.button("Stop Pipeline"):
            st.session_state.pipeline_status = "Stopped"
            st.success("Pipeline stopped successfully.")
        
        return
    
    # Show configuration summary
    st.markdown("### Pipeline Configuration Summary")
    
    # Input sequences
    st.markdown("**Input Sequences:**")
    st.markdown(f"- Hemagglutinin Accession: `{st.session_state.config['ha_accession']}`")
    st.markdown(f"- Neuraminidase Accession: `{st.session_state.config['na_accession']}`")
    
    # Pipeline settings
    st.markdown("**Pipeline Settings:**")
    st.markdown(f"- Output Directory: `{st.session_state.config['outdir']}`")
    st.markdown(f"- Experiment ID: `{st.session_state.config['experiment_id']}`")
    st.markdown(f"- Run Molecular Dynamics: `{'Yes' if st.session_state.config['run_md'] else 'No'}`")
    
    # Calculate experiment output directory
    experiment_output = f"{st.session_state.config['outdir']}/results_{st.session_state.config['experiment_id']}"
    
    # Run the pipeline
    st.markdown("### Execute Pipeline")
    
    # Option to show command preview
    show_command = st.checkbox("Show Command Preview", value=False)
    
    # Construct the command
    cmd = "./run_pipeline.sh"
    cmd += f" --ha_accession {st.session_state.config['ha_accession']}"
    cmd += f" --na_accession {st.session_state.config['na_accession']}"
    cmd += f" --outdir {st.session_state.config['outdir']}"
    cmd += f" --experiment_id {st.session_state.config['experiment_id']}"
    
    if st.session_state.config['run_md']:
        cmd += " --run_md"
    
    # Show command preview if requested
    if show_command:
        st.code(cmd, language="bash")
    
    # Run button
    col1, col2 = st.columns([1, 2])
    
    with col1:
        if st.button("Run Pipeline", key="run_pipeline_button"):
            # For demo purposes, we'll use the demo data instead of actually running
            # the pipeline since we don't have the actual pipeline executable
            
            # Create a queue for communication with the thread
            result_queue = queue.Queue()
            st.session_state.result_queue = result_queue
            
            # Start thread to run pipeline asynchronously
            thread = threading.Thread(
                target=run_pipeline_async,
                args=(cmd, result_queue)
            )
            thread.daemon = True
            thread.start()
            
            st.session_state.pipeline_status = "Running"
            st.success("Pipeline started!")
            
            # Rerun to show progress
            st.rerun()
    
    with col2:
        if st.button("Run Demo (Skip Pipeline)", key="run_demo_button"):
            # Load sample results immediately
            st.session_state.results = load_sample_results()
            st.session_state.pipeline_status = "Completed"
            st.success("Demo results loaded successfully! Go to the View Results page to explore.")

def show_results_page():
    st.markdown('<h2 class="sub-header">View Results</h2>', unsafe_allow_html=True)
    
    # Check if results are available
    if not st.session_state.results or all(len(v) == 0 for v in st.session_state.results.values()):
        st.markdown("""
        <div class="info-box">
        No results available yet. Please run the pipeline first or use the "Run Demo" option.
        </div>
        """, unsafe_allow_html=True)
        
        if st.button("Load Demo Results"):
            st.session_state.results = load_sample_results()
            st.success("Demo results loaded successfully!")
            st.rerun()
        
        return
    
    # Create tabs for different result types
    result_type = st.radio("Result Type", [
        "B-cell Epitopes", 
        "T-cell Class I Epitopes", 
        "T-cell Class II Epitopes", 
        "Combined Epitopes", 
        "Vaccine Constructs"
    ])
    
    # Select protein
    protein_options = {
        "B-cell Epitopes": [item['protein'] for item in st.session_state.results['bcell_epitopes']],
        "T-cell Class I Epitopes": [item['protein'] for item in st.session_state.results['tcell_i_epitopes']],
        "T-cell Class II Epitopes": [item['protein'] for item in st.session_state.results['tcell_ii_epitopes']],
        "Combined Epitopes": [item['protein'] for item in st.session_state.results['combined_epitopes']],
        "Vaccine Constructs": [item['protein'] for item in st.session_state.results['vaccine_constructs']]
    }
    
    current_proteins = protein_options.get(result_type, [])
    
    if not current_proteins:
        st.warning(f"No {result_type} results available.")
        return
    
    selected_protein = st.selectbox("Select Protein", current_proteins)
    
    # Display the selected result type
    if result_type == "B-cell Epitopes":
        show_bcell_results(st.session_state.results, selected_protein)
    elif result_type == "T-cell Class I Epitopes":
        show_tcell_i_results(st.session_state.results, selected_protein)
    elif result_type == "T-cell Class II Epitopes":
        show_tcell_ii_results(st.session_state.results, selected_protein)
    elif result_type == "Combined Epitopes":
        show_combined_results(st.session_state.results, selected_protein)
    elif result_type == "Vaccine Constructs":
        show_vaccine_results(st.session_state.results, selected_protein)

def show_bcell_results(results, selected_protein):
    st.markdown(f"### {selected_protein.capitalize()} B-cell Epitope Results")
    
    # Find the selected protein's data
    data = None
    for item in results['bcell_epitopes']:
        if item['protein'] == selected_protein:
            data = item['data']
            break
    
    if data is None:
        st.warning(f"No B-cell epitope data found for {selected_protein}.")
        return
    
    # Create columns for table and chart
    col1, col2 = st.columns([1, 1])
    
    with col1:
        # Display the data table
        st.subheader("Epitope Table")
        st.dataframe(data)
        
        # Provide download link
        csv = data.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name=f"{selected_protein}_bcell_epitopes.csv",
            mime="text/csv"
        )
    
    with col2:
        # Create a bar chart of epitope scores
        if 'Score' in data.columns and 'Epitope' in data.columns:
            st.subheader("Epitope Scores")
            
            # Sort by score
            sorted_data = data.sort_values('Score', ascending=False).head(10)
            
            fig = px.bar(
                sorted_data,
                x='Epitope',
                y='Score',
                title=f"Top B-cell Epitope Scores",
                color='Score',
                color_continuous_scale='greens'
            )
            
            # Update layout
            fig.update_layout(
                xaxis_title="Epitope",
                yaxis_title="Prediction Score",
                yaxis=dict(range=[0, 1])
            )
            
            st.plotly_chart(fig)
    
    # Show the top epitopes highlighted
    st.subheader("Top B-cell Epitope Candidates")
    
    # Get top 3 epitopes
    top_epitopes = data.sort_values('Score', ascending=False).head(3)
    
    for i, row in top_epitopes.iterrows():
        score = row['Score']
        epitope = row['Epitope']
        position = row['Position']
        
        # Color coding based on score
        if score > 0.98:
            color = "#2e7d32"  # Dark green
        elif score > 0.95:
            color = "#43a047"  # Medium green
        else:
            color = "#81c784"  # Light green
        
        st.markdown(f"""
        <div style="background-color: #e8f5e9; padding: 10px; border-left: 5px solid {color}; margin-bottom: 10px;">
            <h4 style="margin: 0; color: {color};">{epitope}</h4>
            <p>Position: {position} | Score: {score:.4f}</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Add paper comparison
    if selected_protein == "hemagglutinin":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) identified the following B-cell epitopes for hemagglutinin:
        
        - **IGTSTLNQR** (Position 216, Score 0.984)
        - **MVSLVKSDQ** (Position 10, Score 0.964)
        
        These were selected based on their ability to bind to T-cell CD4+ (class II HLA).
        """)
    elif selected_protein == "neuraminidase":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) identified the following B-cell epitope for neuraminidase:
        
        - **YNGIITDTI** (Position 188, Score 0.975)
        
        This was selected based on its ability to bind to T-cell CD4+ (class II HLA).
        """)

def show_tcell_i_results(results, selected_protein):
    st.markdown(f"### {selected_protein.capitalize()} T-cell Class I Epitope Results")
    
    # Find the selected protein's data
    data = None
    for item in results['tcell_i_epitopes']:
        if item['protein'] == selected_protein:
            data = item['data']
            break
    
    if data is None:
        st.warning(f"No T-cell Class I epitope data found for {selected_protein}.")
        return
    
    # Create columns for table and chart
    col1, col2 = st.columns([1, 1])
    
    with col1:
        # Display the data table
        st.subheader("Epitope Table")
        st.dataframe(data)
        
        # Provide download link
        csv = data.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name=f"{selected_protein}_tcell_i_epitopes.csv",
            mime="text/csv"
        )
    
    with col2:
        # Create a grouped bar chart by allele
        if 'Allele' in data.columns:
            st.subheader("Epitopes by HLA Allele")
            
            # Count epitopes per allele
            allele_counts = data['Allele'].value_counts().reset_index()
            allele_counts.columns = ['Allele', 'Count']
            
            fig = px.bar(
                allele_counts,
                x='Allele',
                y='Count',
                title="Epitope Distribution by HLA Allele",
                color='Count',
                color_continuous_scale='greens'
            )
            
            fig.update_layout(
                xaxis_title="HLA Allele",
                yaxis_title="Number of Epitopes"
            )
            
            st.plotly_chart(fig)
    
    # Show the top epitopes highlighted
    st.subheader("Top T-cell Class I Epitope Candidates")
    
    # Get epitopes sorted by affinity and score
    if 'Affinity' in data.columns and 'Score' in data.columns:
        # Filter strong binders with high scores
        strong_binders = data[data['Affinity'] == 'Strong Binder'].sort_values('Score', ascending=False).head(5)
        
        for i, row in strong_binders.iterrows():
            score = row['Score']
            epitope = row['Epitope']
            position = row['Position']
            allele = row['Allele']
            
            st.markdown(f"""
            <div style="background-color: #e8f5e9; padding: 10px; border-left: 5px solid #2e7d32; margin-bottom: 10px;">
                <h4 style="margin: 0; color: #2e7d32;">{epitope}</h4>
                <p>Position: {position} | Allele: {allele} | Score: {score:.4f} | Strong Binder</p>
            </div>
            """, unsafe_allow_html=True)
    
    # Add paper comparison
    if selected_protein == "hemagglutinin":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) identified the following T-cell Class I epitopes for hemagglutinin:
        
        - **MEKIVLLLA** (Position 1, ΔG = -40.0148 kcal/mol, HLA-B*35:01)
        - **CPYLGSPSF** (Position 151, ΔG = -28.0471 kcal/mol, HLA-B*07:02)  
        - **KCQTPMGAI** (Position 293, ΔG = -33.6109 kcal/mol, HLA-B*07:02)
        - **KAVDGVTNK** (Position 389, ΔG = -48.9974 kcal/mol, HLA-A*11:01)
        
        These epitopes were selected based on molecular dynamics and binding energy calculations.
        """)
    elif selected_protein == "neuraminidase":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) identified the following T-cell Class I epitopes for neuraminidase:
        
        - **NPNQKIITI** (Position 2, ΔG = -23.8779 kcal/mol, HLA-B*07:02)
        - **CYPDAGEIT** (Position 261, ΔG = -33.0200 kcal/mol, HLA-A*24)
        
        These epitopes were selected based on molecular dynamics and binding energy calculations.
        """)

def show_tcell_ii_results(results, selected_protein):
    st.markdown(f"### {selected_protein.capitalize()} T-cell Class II Epitope Results")
    
    # Find the selected protein's data
    data = None
    for item in results['tcell_ii_epitopes']:
        if item['protein'] == selected_protein:
            data = item['data']
            break
    
    if data is None:
        st.warning(f"No T-cell Class II epitope data found for {selected_protein}.")
        return
    
    # Create columns for table and chart
    col1, col2 = st.columns([1, 1])
    
    with col1:
        # Display the data table
        st.subheader("Epitope Table")
        st.dataframe(data)
        
        # Provide download link
        csv = data.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name=f"{selected_protein}_tcell_ii_epitopes.csv",
            mime="text/csv"
        )
    
    with col2:
        # Create a grouped bar chart by allele
        if 'Allele' in data.columns:
            st.subheader("Epitopes by HLA Allele")
            
            # Count epitopes per allele
            allele_counts = data['Allele'].value_counts().reset_index()
            allele_counts.columns = ['Allele', 'Count']
            
            fig = px.bar(
                allele_counts,
                x='Allele',
                y='Count',
                title="Epitope Distribution by HLA Allele",
                color='Count',
                color_continuous_scale='greens'
            )
            
            fig.update_layout(
                xaxis_title="HLA Allele",
                yaxis_title="Number of Epitopes"
            )
            
            st.plotly_chart(fig)
    
    # Show the top epitopes highlighted
    st.subheader("Top T-cell Class II Epitope Candidates")
    
    # Get epitopes sorted by affinity and score
    if len(data) > 0:
        # Show all epitopes since there are usually few of them
        for i, row in data.iterrows():
            epitope = row['Epitope']
            position = row['Position']
            allele = row['Allele']
            
            # Check if we have affinity information
            affinity_info = ""
            if 'Affinity' in data.columns:
                affinity = row['Affinity']
                affinity_info = f"| {affinity}"
            
            # Check if we have score information
            score_info = ""
            if 'Score' in data.columns:
                score = row['Score']
                score_info = f"| Score: {score:.4f}"
            
            st.markdown(f"""
            <div style="background-color: #e8f5e9; padding: 10px; border-left: 5px solid #43a047; margin-bottom: 10px;">
                <h4 style="margin: 0; color: #43a047;">{epitope}</h4>
                <p>Position: {position} | Allele: {allele} {score_info} {affinity_info}</p>
            </div>
            """, unsafe_allow_html=True)
    
    # Add paper comparison
    if selected_protein == "hemagglutinin":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) identified the following T-cell Class II epitopes for hemagglutinin:
        
        - **MVSLVKSDQ** (Position 10, ΔG = -11.7756 kcal/mol, HLA-DRB1*03:01)
        - **IGTSTLNQR** (Position 216, ΔG = -56.9580 kcal/mol, HLA-DRB1*03:01)
        
        These epitopes were selected based on molecular dynamics and binding energy calculations.
        """)
    elif selected_protein == "neuraminidase":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) identified the following T-cell Class II epitope for neuraminidase:
        
        - **YNGIITDTI** (Position 188, ΔG = -17.8049 kcal/mol, HLA-DRB1*01:01)
        
        This epitope was selected based on molecular dynamics and binding energy calculations.
        """)

def show_combined_results(results, selected_protein):
    st.markdown(f"### {selected_protein.capitalize()} Combined Epitope Results")
    
    # Find the selected protein's data
    data = None
    for item in results['combined_epitopes']:
        if item['protein'] == selected_protein:
            data = item['data']
            break
    
    if data is None:
        st.warning(f"No combined epitope data found for {selected_protein}.")
        return
    
    # Create columns for table and chart
    col1, col2 = st.columns([1, 1])
    
    with col1:
        # Display the data table
        st.subheader("Combined Epitopes")
        st.dataframe(data)
        
        # Provide download link
        csv = data.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name=f"{selected_protein}_combined_epitopes.csv",
            mime="text/csv"
        )
    
    with col2:
        # Create a pie chart of epitope types
        if 'Type' in data.columns:
            st.subheader("Epitope Type Distribution")
            
            # Count epitopes by type
            type_counts = data['Type'].value_counts().reset_index()
            type_counts.columns = ['Type', 'Count']
            
            fig = px.pie(
                type_counts,
                names='Type',
                values='Count',
                title="Epitope Types",
                color_discrete_sequence=px.colors.sequential.Greens
            )
            
            st.plotly_chart(fig)
    
    # Binding energy visualization if available
    if 'ΔG binding' in data.columns:
        st.subheader("Binding Energy Comparison")
        
        # Sort by binding energy
        sorted_data = data.sort_values('ΔG binding')
        
        fig = px.bar(
            sorted_data,
            x='Epitope',
            y='ΔG binding',
            title="Epitope Binding Energy (ΔG)",
            color='ΔG binding',
            color_continuous_scale='greens_r'  # Reversed to make stronger bindings darker
        )
        
        fig.update_layout(
            xaxis_title="Epitope",
            yaxis_title="Binding Energy (kcal/mol)"
        )
        
        st.plotly_chart(fig)
    
    # Show the final selected epitopes
    st.subheader("Selected Epitopes for Vaccine Design")
    
    if 'ΔG binding' in data.columns:
        # Sort by binding energy (stronger binding first)
        selected_epitopes = data.sort_values('ΔG binding').head(4)
    else:
        # Just show all epitopes if no binding energy data
        selected_epitopes = data
    
    for i, row in selected_epitopes.iterrows():
        epitope = row['Epitope']
        position = row['Position']
        epitope_type = row['Type'] if 'Type' in data.columns else "Unknown"
        
        # Get binding energy if available
        binding_info = ""
        if 'ΔG binding' in data.columns:
            binding = row['ΔG binding']
            binding_info = f"| ΔG: {binding:.4f} kcal/mol"
        
        # Get HLA info if available
        hla_info = ""
        if 'HLA' in data.columns:
            hla = row['HLA']
            hla_info = f"| HLA: {hla}"
        
        st.markdown(f"""
        <div style="background-color: #e8f5e9; padding: 10px; border-left: 5px solid #1b5e20; margin-bottom: 10px;">
            <h4 style="margin: 0; color: #1b5e20;">{epitope}</h4>
            <p>Position: {position} | Type: {epitope_type} {hla_info} {binding_info}</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Add paper comparison
    if selected_protein == "hemagglutinin":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) selected the following epitopes for hemagglutinin vaccine design:
        
        - **MEKIVLLLA** (ΔG = -40.0148 kcal/mol, T-cell I)
        - **CPYLGSPSF** (ΔG = -28.0471 kcal/mol, T-cell I)
        - **KCQTPMGAI** (ΔG = -33.6109 kcal/mol, T-cell I)
        - **IGTSTLNQR** (ΔG = -56.9580 kcal/mol, B-cell & T-cell II)
        
        These epitopes were selected based on binding energy calculations and molecular dynamics simulations.
        """)
    elif selected_protein == "neuraminidase":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) selected the following epitopes for neuraminidase vaccine design:
        
        - **NPNQKIITI** (ΔG = -23.8779 kcal/mol, T-cell I)
        - **CYPDAGEIT** (ΔG = -33.0200 kcal/mol, T-cell I)
        
        These epitopes were selected based on binding energy calculations and molecular dynamics simulations.
        """)

def show_vaccine_results(results, selected_protein):
    st.markdown(f"### {selected_protein.capitalize()} Vaccine Construct")
    
    # Find the selected protein's data
    content = None
    for item in results['vaccine_constructs']:
        if item['protein'] == selected_protein:
            content = item['content']
            break
    
    if content is None:
        st.warning(f"No vaccine construct found for {selected_protein}.")
        return
    
    # Display the FASTA content
    st.subheader("Vaccine Sequence (FASTA)")
    
    # Format FASTA for display
    st.markdown(f'<div class="terminal">{content}</div>', unsafe_allow_html=True)
    
    # Extract sequence
    sequence = ""
    for line in content.split('\n'):
        if not line.startswith('>'):
            sequence += line.strip()
    
    # Create columns for properties and visualization
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("Sequence Properties")
        
        # Sequence length
        st.metric("Sequence Length", len(sequence))
        
        # Amino acid composition
        aa_counts = {}
        for aa in sorted(set(sequence)):
            aa_counts[aa] = sequence.count(aa)
        
        aa_df = pd.DataFrame({
            'Amino Acid': list(aa_counts.keys()),
            'Count': list(aa_counts.values())
        })
        
        fig = px.bar(
            aa_df,
            x='Amino Acid',
            y='Count',
            title='Amino Acid Composition',
            color='Count',
            color_continuous_scale='greens'
        )
        
        st.plotly_chart(fig)
    
    with col2:
        st.subheader("Epitope Visualization")
        
        # Detect GPGPG linkers
        epitopes = sequence.split("GPGPG")
        
        # Create a visualization
        fig = go.Figure()
        
        y_pos = 1
        current_pos = 0
        
        for i, epitope in enumerate(epitopes):
            if epitope:
                # Add epitope segment
                epitope_len = len(epitope)
                fig.add_trace(go.Scatter(
                    x=[current_pos, current_pos + epitope_len],
                    y=[y_pos, y_pos],
                    mode='lines',
                    line=dict(width=20, color='#2e7d32'),
                    name=f'Epitope {i+1}'
                ))
                
                # Add label
                fig.add_annotation(
                    x=current_pos + epitope_len/2,
                    y=y_pos + 0.2,
                    text=f'E{i+1}',
                    showarrow=False,
                    font=dict(color='white', size=12)
                )
                
                current_pos += epitope_len
                
                # Add linker if not the last epitope
                if i < len(epitopes) - 1:
                    linker_len = 5  # GPGPG
                    fig.add_trace(go.Scatter(
                        x=[current_pos, current_pos + linker_len],
                        y=[y_pos, y_pos],
                        mode='lines',
                        line=dict(width=20, color='#ffa000'),
                        name='Linker'
                    ))
                    
                    # Add label
                    fig.add_annotation(
                        x=current_pos + linker_len/2,
                        y=y_pos + 0.2,
                        text='L',
                        showarrow=False,
                        font=dict(color='white', size=12)
                    )
                    
                    current_pos += linker_len
        
        # Update layout
        fig.update_layout(
            xaxis=dict(title='Position', showgrid=False),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
            showlegend=False,
            margin=dict(l=20, r=20, t=20, b=20),
            height=200
        )
        
        st.plotly_chart(fig)
        
        # Show epitope details
        st.subheader("Epitope Details")
        
        for i, epitope in enumerate(epitopes):
            if epitope:
                st.markdown(f"**Epitope {i+1}:** `{epitope}` (Length: {len(epitope)})")
    
    # Provide download link
    st.download_button(
        label=f"Download {selected_protein} vaccine construct",
        data=content,
        file_name=f"{selected_protein}_vaccine_construct.fasta",
        mime="text/plain"
    )
    
    # Show paper comparison
    st.subheader("Notes on Vaccine Design")
    
    if selected_protein == "hemagglutinin":
        st.markdown("""
        The vaccine construct uses GPGPG linkers between epitopes to prevent junctional epitopes and enhance presentation.
        
        In the paper by Tambunan et al. (2016), the hemagglutinin epitopes were selected based on:
        1. Strong binding affinity to HLA molecules
        2. Conservation across H5N1 strains
        3. Stable molecular dynamics simulation results
        
        The selected epitopes (MEKIVLLLA, CPYLGSPSF, KCQTPMGAI, IGTSTLNQR) showed binding energies 
        from -28.05 to -56.96 kcal/mol, indicating strong and stable binding to their respective HLA molecules.
        """)
    elif selected_protein == "neuraminidase":
        st.markdown("""
        The vaccine construct uses GPGPG linkers between epitopes to prevent junctional epitopes and enhance presentation.
        
        In the paper by Tambunan et al. (2016), the neuraminidase epitopes were selected based on:
        1. Strong binding affinity to HLA molecules
        2. Conservation across H5N1 strains
        3. Stable molecular dynamics simulation results
        
        The selected epitopes (NPNQKIITI, CYPDAGEIT) showed binding energies 
        from -23.88 to -33.02 kcal/mol, indicating strong and stable binding to their respective HLA molecules.
        """)
    elif selected_protein == "combined":
        st.markdown("""
        This combined vaccine construct incorporates epitopes from both hemagglutinin and neuraminidase proteins, 
        providing broader coverage against the H5N1 virus. The GPGPG linkers separate epitopes to prevent junctional 
        epitopes and enhance presentation.
        
        By targeting both major surface proteins of the H5N1 virus, this multi-epitope vaccine design aims to:
        1. Stimulate both B-cell and T-cell responses
        2. Provide broader coverage against viral variants
        3. Improve protection compared to single-protein vaccines
        
        This approach follows recommendations from Tambunan et al. (2016) for developing effective epitope-based vaccines.
        """)

def show_about_page():
    st.markdown('<h2 class="sub-header">About</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
    <h3>H5N1 Epitope-Based Vaccine Design Application</h3>
    <p>This application implements the methodology described in Tambunan et al. (2016) for designing epitope-based vaccines against H5N1 avian influenza virus.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Paper details
    st.subheader("Research Paper")
    
    st.markdown("""
    **Title:** Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions
    
    **Authors:** Usman Sumo Friend Tambunan, Feimmy Ruth Pratiwi Sipahutar, Arli Aditya Parikesit, and Djati Kerami
    
    **Journal:** Bioinformatics and Biology Insights 2016:10 27-35
    
    **DOI:** 10.4137/BBI.S38378
    
    [View Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4841451/)
    """)
    
    # Methodology
    st.subheader("Methodology")
    
    st.markdown("""
    The methodology implemented in this application follows these key steps:
    
    1. **Sequence Retrieval**: Hemagglutinin and neuraminidase protein sequences from H5N1 Indonesian strain are used as input.
    
    2. **Epitope Prediction**:
       - B-cell epitopes are predicted using BCPREDS with 80% specificity and 9-mer length
       - T-cell Class I epitopes are predicted using ProPred I and NetMHCpan
       - T-cell Class II epitopes are predicted using ProPred and NetMHCIIpan
    
    3. **Epitope Selection**:
       - B-cell epitopes that can bind to T-cell CD4+ are identified
       - T-cell epitopes are filtered based on proteasomal cleavage and TAP binding
       - HLA binding affinity is evaluated using multiple prediction methods
    
    4. **Binding Energy Calculation**:
       - Molecular dynamics simulation to calculate binding energy (ΔG) between epitopes and HLA molecules
       - Epitopes with strongest binding (most negative ΔG values) are selected
    
    5. **Vaccine Design**:
       - Selected epitopes are combined with GPGPG linkers to create multi-epitope vaccine constructs
    """)
    
    # Key findings
    st.subheader("Key Findings")
    
    st.markdown("""
    The key findings from the paper that guided this application:
    
    - **Hemagglutinin Epitopes**: Four epitopes (MEKIVLLLA, CPYLGSPSF, KCQTPMGAI, IGTSTLNQR) were identified as having the best binding affinity to HLA molecules.
    
    - **Neuraminidase Epitopes**: Two epitopes (NPNQKIITI, CYPDAGEIT) were identified as having the best binding affinity to HLA molecules.
    
    - **Binding Energies**: The strongest binding was observed for IGTSTLNQR with HLA-DRB1*0301 (ΔG = -56.958 kcal/mol).
    
    - **Molecular Dynamics**: The epitope-HLA structures stabilized after 100 ps of simulation, confirming their structural stability.
    """)
    
    # Application information
    st.subheader("Application Information")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        **Version:** 1.0.0
        
        **Released:** April 2025
        
        **Developed by:** Kalbe Bioinformatics Team
        
        **Based on:** Tambunan et al. (2016)
        """)
    
    with col2:
        st.markdown("""
        **Technologies:**
        - NextFlow (Pipeline orchestration)
        - Streamlit (User interface)
        - Bioinformatics tools (BCPREDS, NetMHCpan, etc.)
        - Molecular dynamics simulation
        """)
    
    # Disclaimer
    st.markdown("""
    <div class="warning-box">
    <h4>Disclaimer</h4>
    <p>This application is provided for research and educational purposes only. It is not intended for clinical use or vaccine development without proper validation through wet-lab experiments and clinical trials.</p>
    </div>
    """, unsafe_allow_html=True)

# Run the main function
if __name__ == "__main__":
    main()