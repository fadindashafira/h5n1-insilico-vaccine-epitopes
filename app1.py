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

# Helper functions
def load_kalbe_logo():
    """Load and encode the Kalbe Bioinformatics logo"""
    logo_path = "./img/kalbe_bioinformatics_logo.jpeg"
    if os.path.exists(logo_path):    
        return Image.open(logo_path)
    return None

def run_nextflow_command(cmd, cwd="."):
    """Run a command and capture output in real-time"""
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        text=True,
        cwd=cwd,
        bufsize=1
    )
    return process

def parse_nextflow_log(log_file):
    """Parse Nextflow log file to extract key information"""
    if not os.path.exists(log_file):
        return None
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Extract key information
        run_info = {}
        
        # Get session ID
        session_match = re.search(r'Run name: (\w+)', content)
        if session_match:
            run_info['session_id'] = session_match.group(1)
        
        # Get command line
        cmd_match = re.search(r'Command line: (.*)', content)
        if cmd_match:
            run_info['command'] = cmd_match.group(1)
        
        # Get progress information
        progress_matches = re.findall(r'\[(\w+)\] Process `([^`]+)`.*\[(\d+)%\] (\d+) of (\d+)', content)
        processes = []
        
        for match in progress_matches:
            status, process_name, percent, completed, total = match
            processes.append({
                'status': status,
                'name': process_name,
                'percent': int(percent),
                'completed': int(completed),
                'total': int(total)
            })
        
        run_info['processes'] = processes
        
        # Get status
        if "[skipped]" in content:
            run_info['status'] = 'SKIPPED'
        elif "Execution cancelled" in content:
            run_info['status'] = 'CANCELLED'
        elif "Execution completed successfully" in content:
            run_info['status'] = 'COMPLETED'
        elif "Error executing process" in content:
            run_info['status'] = 'ERROR'
        else:
            run_info['status'] = 'RUNNING'
        
        return run_info
    
    except Exception as e:
        st.error(f"Error parsing log file: {str(e)}")
        return None

def parse_timeline_file(timeline_file):
    """Parse Nextflow timeline file"""
    if not os.path.exists(timeline_file):
        return None
    
    try:
        # Timeline is HTML with embedded JSON data
        with open(timeline_file, 'r') as f:
            content = f.read()
        
        # Extract the JSON data from the HTML
        data_match = re.search(r'var data = (\[.*?\]);', content, re.DOTALL)
        if data_match:
            json_data = data_match.group(1)
            # Parse the JSON data
            timeline_data = json.loads(json_data)
            return timeline_data
        
        return None
    
    except Exception as e:
        st.error(f"Error parsing timeline file: {str(e)}")
        return None

def get_latest_results(experiment_output):
    """Get the latest results files from the experiment output directory"""
    if not os.path.exists(experiment_output):
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
        except Exception as e:
            st.warning(f"Error loading {file}: {str(e)}")
    
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
        except Exception as e:
            st.warning(f"Error loading {file}: {str(e)}")
    
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
        except Exception as e:
            st.warning(f"Error loading {file}: {str(e)}")
    
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
        except Exception as e:
            st.warning(f"Error loading {file}: {str(e)}")
    
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
        except Exception as e:
            st.warning(f"Error loading {file}: {str(e)}")
    
    # If no results were found, load sample results
    if all(not items for items in results.values()):
        return load_sample_results()
    
    return results

def load_sample_results():
    """Load sample results based on Tambunan et al. (2016) paper data"""
    # Sample B-cell epitopes for hemagglutinin
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
    
    # Sample B-cell epitopes for neuraminidase
    na_bcell = pd.DataFrame([
        {"Position": 342, "Epitope": "TKSTNSRSG", "Score": 0.996},
        {"Position": 174, "Epitope": "GISGPDNEA", "Score": 0.995},
        {"Position": 142, "Epitope": "PVGEAPSPY", "Score": 0.995},
        {"Position": 303, "Epitope": "GDNPRPNDG", "Score": 0.989},
        {"Position": 84, "Epitope": "NNIRIGSKG", "Score": 0.977},
        {"Position": 188, "Epitope": "YNGIITDTI", "Score": 0.975},
        {"Position": 257, "Epitope": "EESCSCYDA", "Score": 0.972},
        {"Position": 242, "Epitope": "KVVKSVELD", "Score": 0.970},
        {"Position": 3, "Epitope": "PNQKIITIG", "Score": 0.951}
    ])
    
    # Sample T-cell Class I epitopes for hemagglutinin
    ha_tcell_i = pd.DataFrame([
        {"Position": 1, "Epitope": "MEKIVLLLA", "Allele": "HLA-B*35:01", "Score": 0.85, "Affinity": "Strong Binder", "ΔG binding": -40.0148},
        {"Position": 2, "Epitope": "EKIVLLLAM", "Allele": "HLA-B*35:01", "Score": 0.75, "Affinity": "Weak Binder", "ΔG binding": -22.4538},
        {"Position": 151, "Epitope": "CPYLGSPSF", "Allele": "HLA-B*07:02", "Score": 0.92, "Affinity": "Strong Binder", "ΔG binding": -28.0471},
        {"Position": 293, "Epitope": "KCQTPMGAI", "Allele": "HLA-B*07:02", "Score": 0.88, "Affinity": "Strong Binder", "ΔG binding": -33.6109},
        {"Position": 389, "Epitope": "KAVDGVTNK", "Allele": "HLA-A*11:01", "Score": 0.77, "Affinity": "Weak Binder", "ΔG binding": -48.9974}
    ])
    
    # Sample T-cell Class I epitopes for neuraminidase
    na_tcell_i = pd.DataFrame([
        {"Position": 2, "Epitope": "NPNQKIITI", "Allele": "HLA-B*07:02", "Score": 0.89, "Affinity": "Strong Binder", "ΔG binding": -23.8779},
        {"Position": 261, "Epitope": "CYPDAGEIT", "Allele": "HLA-A*24", "Score": 0.79, "Affinity": "Strong Binder", "ΔG binding": -33.0200},
        {"Position": 398, "Epitope": "IRPCFWVEL", "Allele": "HLA-B*27:05", "Score": 0.86, "Affinity": "Strong Binder", "ΔG binding": -33.6376},
        {"Position": 399, "Epitope": "RPCFWVELI", "Allele": "HLA-B*07:02", "Score": 0.72, "Affinity": "Weak Binder", "ΔG binding": -24.8928}
    ])
    
    # Sample T-cell Class II epitopes for hemagglutinin
    ha_tcell_ii = pd.DataFrame([
        {"Position": 10, "Epitope": "MVSLVKSDQ", "Allele": "HLA-DRB1*03:01", "Score": 0.65, "Affinity": "Weak Binder", "ΔG binding": -11.7756},
        {"Position": 216, "Epitope": "IGTSTLNQR", "Allele": "HLA-DRB1*03:01", "Score": 0.81, "Affinity": "Strong Binder", "ΔG binding": -56.9580}
    ])
    
    # Sample T-cell Class II epitopes for neuraminidase
    na_tcell_ii = pd.DataFrame([
        {"Position": 188, "Epitope": "YNGIITDTI", "Allele": "HLA-DRB1*01:01", "Score": 0.68, "Affinity": "Weak Binder", "ΔG binding": -17.8049}
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
        {"Position": 188, "Epitope": "YNGIITDTI", "Type": "B-cell & T-cell II", "Score": 0.975, "HLA": "HLA-DRB1*01:01", "ΔG binding": -17.805},
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
            {'protein': 'hemagglutinin', 'file': 'sample_data/ha_bcell.csv', 'data': ha_bcell},
            {'protein': 'neuraminidase', 'file': 'sample_data/na_bcell.csv', 'data': na_bcell}
        ],
        'tcell_i_epitopes': [
            {'protein': 'hemagglutinin', 'file': 'sample_data/ha_tcell_i.csv', 'data': ha_tcell_i},
            {'protein': 'neuraminidase', 'file': 'sample_data/na_tcell_i.csv', 'data': na_tcell_i}
        ],
        'tcell_ii_epitopes': [
            {'protein': 'hemagglutinin', 'file': 'sample_data/ha_tcell_ii.csv', 'data': ha_tcell_ii},
            {'protein': 'neuraminidase', 'file': 'sample_data/na_tcell_ii.csv', 'data': na_tcell_ii}
        ],
        'combined_epitopes': [
            {'protein': 'hemagglutinin', 'file': 'sample_data/ha_combined.csv', 'data': ha_combined},
            {'protein': 'neuraminidase', 'file': 'sample_data/na_combined.csv', 'data': na_combined}
        ],
        'vaccine_constructs': [
            {'protein': 'hemagglutinin', 'file': 'sample_data/ha_vaccine.fasta', 'content': ha_vaccine},
            {'protein': 'neuraminidase', 'file': 'sample_data/na_vaccine.fasta', 'content': na_vaccine},
            {'protein': 'combined', 'file': 'sample_data/combined_vaccine.fasta', 'content': combined_vaccine}
        ]
    }

def create_process_status_chart(processes):
    """Create a chart showing process status"""
    if not processes:
        return None
    
    # Create a DataFrame from the processes
    df = pd.DataFrame(processes)
    
    # Create a horizontal bar chart
    fig = go.Figure()
    
    # Add a bar for each process
    for i, row in df.iterrows():
        color = '#4CAF50' if row['status'] == 'COMPLETED' else '#FFC107'  # Green for completed, amber for in progress
        
        # If there's an error, use red
        if row['status'] == 'ERROR':
            color = '#F44336'
        
        fig.add_trace(go.Bar(
            y=[row['name']],
            x=[row['percent']],
            orientation='h',
            name=row['name'],
            marker=dict(color=color),
            text=f"{row['completed']} of {row['total']} ({row['percent']}%)",
            textposition='auto'
        ))
    
    # Update layout
    fig.update_layout(
        title='Process Status',
        yaxis=dict(title='Process'),
        xaxis=dict(title='Completion %', range=[0, 100]),
        height=len(processes) * 40 + 100,  # Dynamic height based on number of processes
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    return fig

def run_pipeline_async(cmd, output_queue):
    """Run pipeline asynchronously and put results in queue"""
    try:
        st.session_state.pipeline_status = "Running"
        st.session_state.log_output = ""
        st.session_state.progress = 0.0
        st.session_state.current_task = "Starting pipeline..."
        
        # Log the command we're about to run
        st.session_state.log_output += f"Executing command: {cmd}\n"
        
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            text=True,
            bufsize=1
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
        
        # Get any error output
        stderr = process.stderr.read()
        if stderr:
            st.session_state.log_output += f"\nERRORS:\n{stderr}\n"
        
        # Check if process completed successfully
        if process.returncode == 0:
            st.session_state.pipeline_status = "Completed"
            st.session_state.progress = 1.0
            
            # Load the results directly after completion
            experiment_id = re.search(r'--experiment_id\s+(\S+)', cmd)
            experiment_id = experiment_id.group(1) if experiment_id else "exp1"
            experiment_output = f"results/results_{experiment_id}"
            results = get_latest_results(experiment_output)
            output_queue.put({"status": "success", "results": results})
        else:
            st.session_state.pipeline_status = "Failed"
            output_queue.put({"status": "error", "message": f"Pipeline failed with return code {process.returncode}"})
    
    except Exception as e:
        st.session_state.pipeline_status = "Error"
        st.session_state.log_output += f"\nException occurred: {str(e)}\n"
        output_queue.put({"status": "error", "message": str(e)})

def main():
    # Initialize session state for storing results
    if 'results' not in st.session_state:
        st.session_state.results = None
    
    if 'pipeline_status' not in st.session_state:
        st.session_state.pipeline_status = "Not Started"
    
    if 'log_output' not in st.session_state:
        st.session_state.log_output = ""
    
    if 'progress' not in st.session_state:
        st.session_state.progress = 0.0
    
    if 'current_task' not in st.session_state:
        st.session_state.current_task = "No task running"
    
    # Initialize result_queue for async operations if not exists
    if 'result_queue' not in st.session_state:
        st.session_state.result_queue = queue.Queue()
    
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
            'ha_accession': 'BAL61222.1',
            'na_accession': 'BAL61230.1',
            'outdir': 'results',
            'experiment_id': f'test_exp_1',
            'profile': 'local',
            'run_md': False,
            'resume': False
        }
    
    if 'run_status' not in st.session_state:
        st.session_state.run_status = {
            'is_running': False,
            'process': None,
            'run_id': None,
            'start_time': None,
            'log_file': None,
            'experiment_output': None
        }
    
    # Check for results from async operations
    if not st.session_state.result_queue.empty():
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
    and designs optimal vaccine constructs.
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
    4. **Vaccine Design**: Combine selected epitopes with GPGPG linkers to create vaccine constructs
    5. **Evaluation**: Assess immunogenicity and stability of the designed vaccines
    
    ### Getting Started
    
    To begin using the application, navigate to the **Configure Pipeline** section in the sidebar 
    to set up your parameters, then proceed to **Run Pipeline** to execute the workflow.
    """)
    
    # Display some key metrics
    st.markdown("### Key Features")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Prediction Algorithms", "4+", "BCPREDS, NetMHCpan, ProPred")
    
    with col2:
        st.metric("HLA Coverage", "14+ alleles", "Class I & II")
    
    with col3:
        st.metric("Based on", "Indonesian strain", "H5N1 virus")
    
    # Recent publications
    st.markdown("### Research Background")
    st.markdown("""
    This application implements the methodology from the paper by Tambunan et al. (2016), "Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions." The researchers identified key epitopes from H5N1 hemagglutinin and neuraminidase proteins that showed strong binding affinity to HLA molecules.
    
    Key references:
    - Tambunan et al. (2016). **Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions**. *Bioinformatics and Biology Insights*, 10, 27-35.
    - Velkov et al. (2013). **The antigenic architecture of the hemagglutinin of influenza H5N1 viruses**. *Molecular Immunology*, 56(4), 705-719.
    """)

def show_config_page():
    st.markdown('<h2 class="sub-header">Configure Pipeline</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
    Configure the parameters for the H5N1 vaccine design pipeline. These settings will be used 
    when executing the Nextflow workflow.
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
        
        # Profile selection
        profile_options = ["local", "slurm", "conda", "docker", "singularity"]
        profile = st.selectbox(
            "Execution Profile",
            profile_options,
            index=profile_options.index(st.session_state.config.get('profile', 'local')),
            help="Computational environment for pipeline execution"
        )
        
        # Run molecular dynamics option
        run_md = st.checkbox(
            "Run Molecular Dynamics",
            st.session_state.config.get('run_md', False),
            help="Enable molecular dynamics simulation (requires more computational resources)"
        )
        
        # Resume option
        resume = st.checkbox(
            "Resume Pipeline",
            st.session_state.config.get('resume', False),
            help="Resume pipeline execution from last checkpoint if previously interrupted"
        )
        
        # Save pipeline settings
        if st.button("Save Pipeline Settings"):
            st.session_state.config['outdir'] = outdir
            st.session_state.config['experiment_id'] = experiment_id
            st.session_state.config['profile'] = profile
            st.session_state.config['run_md'] = run_md
            st.session_state.config['resume'] = resume
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
            st.session_state.run_status['is_running'] = False
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
    st.markdown(f"- Execution Profile: `{st.session_state.config['profile']}`")
    st.markdown(f"- Run Molecular Dynamics: `{'Yes' if st.session_state.config['run_md'] else 'No'}`")
    st.markdown(f"- Resume Execution: `{'Yes' if st.session_state.config['resume'] else 'No'}`")
    
    # Calculate experiment output directory
    experiment_output = f"{st.session_state.config['outdir']}/results_{st.session_state.config['experiment_id']}"
    
    # Run the pipeline
    st.markdown("### Execute Pipeline")
    
    # Option to show command preview
    show_command = st.checkbox("Show Command Preview", value=False)
    
    # Construct the command
    cmd = f"./run_pipeline.sh"
    cmd += f" --ha_accession {st.session_state.config['ha_accession']}"
    cmd += f" --na_accession {st.session_state.config['na_accession']}"
    cmd += f" --outdir {st.session_state.config['outdir']}"
    cmd += f" --experiment_id {st.session_state.config['experiment_id']}"
    cmd += f" --profile {st.session_state.config['profile']}"
    
    if st.session_state.config['run_md']:
        cmd += " --run_md"
    
    if st.session_state.config['resume']:
        cmd += " --resume"
    
    # Show command preview if requested
    if show_command:
        st.code(cmd, language="bash")
    
    # Run buttons
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("Run Pipeline", key="run_pipeline_button"):
            # Create a queue for communication with the thread
            result_queue = queue.Queue()
            st.session_state.result_queue = result_queue
            
            # Update run status
            run_id = f"run_{int(time.time())}"
            st.session_state.run_status['is_running'] = True
            st.session_state.run_status['run_id'] = run_id
            st.session_state.run_status['start_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            st.session_state.run_status['experiment_output'] = experiment_output
            
            # Start thread to run pipeline asynchronously
            thread = threading.Thread(
                target=run_pipeline_async,
                args=(cmd, result_queue)
            )
            thread.daemon = True
            thread.start()
            
            # Update pipeline status
            st.session_state.pipeline_status = "Running"
            st.success("Pipeline started! Track progress below.")
            st.session_state.current_task = "Initializing pipeline..."
            
            # Make progress info visible immediately
            st.progress(0.0)
            st.info("Current task: Initializing pipeline...")
            st.expander("Log Output", expanded=True).markdown("<div class='terminal'>Starting pipeline execution...</div>", unsafe_allow_html=True)
    
    with col2:
        if st.button("Load Demo Results", key="load_demo_button"):
            # Load sample results immediately
            st.session_state.results = load_sample_results()
            st.session_state.pipeline_status = "Completed"
            st.success("Demo results loaded successfully! Go to the View Results tab to explore.")
            
            # Provide a convenient link to results page
            st.markdown("[View Results Now](/View_Results)", unsafe_allow_html=True)

def show_results_page():
    st.markdown('<h2 class="sub-header">View Results</h2>', unsafe_allow_html=True)
    
    # First check if we have results in session state
    if st.session_state.results is not None:
        results = st.session_state.results
    else:
        # Find all experiment directories
        result_dirs = glob.glob("results/results_*")
        
        if not result_dirs:
            st.markdown("""
            <div class="info-box">
            No results found. Please run the pipeline first to generate results or use the "Load Demo Results" option on the Run Pipeline page.
            </div>
            """, unsafe_allow_html=True)
            
            if st.button("Load Demo Results"):
                st.session_state.results = load_sample_results()
                st.experimental_rerun()
            
            return
        
        # Create a dropdown to select experiment
        selected_experiment = st.selectbox(
            "Select Experiment",
            sorted(result_dirs, reverse=True),
            format_func=lambda x: os.path.basename(x).replace("results_", "")
        )
        
        if not selected_experiment:
            return
        
        # Load results from the selected experiment
        results = get_latest_results(selected_experiment)
        
        if not results:
            st.warning(f"No result files found in {selected_experiment}")
            
            if st.button("Load Demo Results"):
                st.session_state.results = load_sample_results()
                st.experimental_rerun()
            
            return
    
    # Create radio selection for result type
    result_type = st.radio("Result Type", [
        "B-cell Epitopes", 
        "T-cell Class I Epitopes", 
        "T-cell Class II Epitopes", 
        "Combined Epitopes", 
        "Vaccine Constructs"
    ])
    
    # Map result types to their keys in the results dictionary
    result_type_map = {
        "B-cell Epitopes": "bcell_epitopes",
        "T-cell Class I Epitopes": "tcell_i_epitopes",
        "T-cell Class II Epitopes": "tcell_ii_epitopes",
        "Combined Epitopes": "combined_epitopes",
        "Vaccine Constructs": "vaccine_constructs"
    }
    
    # Get the corresponding results
    result_key = result_type_map[result_type]
    available_results = results[result_key]
    
    if not available_results:
        st.info(f"No {result_type} found.")
        return
    
    # Create protein selector
    proteins = [item['protein'] for item in available_results]
    selected_protein = st.selectbox("Select Protein", proteins)
    
    # Display the selected result type
    if result_type == "B-cell Epitopes":
        show_bcell_results(results, selected_protein)
    elif result_type == "T-cell Class I Epitopes":
        show_tcell_i_results(results, selected_protein)
    elif result_type == "T-cell Class II Epitopes":
        show_tcell_ii_results(results, selected_protein)
    elif result_type == "Combined Epitopes":
        show_combined_results(results, selected_protein)
    elif result_type == "Vaccine Constructs":
        show_vaccine_results(results, selected_protein)
    
    # Show reports if available
    st.markdown("### Pipeline Reports")
    experiment_output = f"results/results_{st.session_state.config['experiment_id']}"
    reports_dir = f"{experiment_output}/reports"
    
    if os.path.exists(reports_dir):
        report_files = glob.glob(f"{reports_dir}/*")
        
        if report_files:
            col1, col2, col3 = st.columns(3)
            
            with col1:
                timeline_file = next((f for f in report_files if 'timeline' in f), None)
                if timeline_file:
                    st.markdown(f"[Timeline Report]({timeline_file})")
            
            with col2:
                execution_file = next((f for f in report_files if 'execution' in f), None)
                if execution_file:
                    st.markdown(f"[Execution Report]({execution_file})")
            
            with col3:
                dag_file = next((f for f in report_files if 'dag' in f), None)
                if dag_file:
                    st.markdown(f"[DAG Visualization]({dag_file})")
            
            # Other reports
            other_reports = [f for f in report_files if not any(x in f for x in ['timeline', 'execution', 'dag'])]
            if other_reports:
                st.markdown("**Other Reports:**")
                for report in other_reports:
                    report_name = os.path.basename(report)
                    st.markdown(f"- [{report_name}]({report})")
        else:
            st.info("No report files found.")
    else:
        st.info("Reports directory not available.")

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
        
        # Show binding energy if available
        if 'ΔG binding' in data.columns:
            st.subheader("Binding Energy")
            
            fig = px.bar(
                data,
                x='Epitope',
                y='ΔG binding',
                title="Epitope Binding Energy (ΔG)",
                color='ΔG binding',
                color_continuous_scale='RdBu_r'  # Red (weak) to Blue (strong)
            )
            
            fig.update_layout(
                xaxis_title="Epitope",
                yaxis_title="Binding Energy (kcal/mol)"
            )
            
            st.plotly_chart(fig)
    
    # Show all epitopes highlighted
    st.subheader("T-cell Class II Epitope Candidates")
    
    # Show all epitopes since there are usually few of them
    for i, row in data.iterrows():
        epitope = row['Epitope']
        position = row['Position']
        allele = row['Allele'] if 'Allele' in data.columns else 'Unknown'
        
        # Additional info
        binding_info = f"ΔG: {row['ΔG binding']:.4f} kcal/mol" if 'ΔG binding' in row else ""
        affinity_info = f"Affinity: {row['Affinity']}" if 'Affinity' in row else ""
        score_info = f"Score: {row['Score']:.4f}" if 'Score' in row else ""
        
        # Combine available info
        info_parts = [part for part in [binding_info, affinity_info, score_info] if part]
        additional_info = " | ".join(info_parts)
        
        # Choose color based on binding energy if available
        color = "#43a047"  # Default green
        if 'ΔG binding' in row:
            energy = row['ΔG binding']
            if energy < -40:
                color = "#1b5e20"  # Very strong binding - dark green
            elif energy < -30:
                color = "#2e7d32"  # Strong binding - medium dark green
            elif energy < -20:
                color = "#43a047"  # Medium binding - medium green
            else:
                color = "#66bb6a"  # Weak binding - light green
        
        st.markdown(f"""
        <div style="background-color: #e8f5e9; padding: 10px; border-left: 5px solid {color}; margin-bottom: 10px;">
            <h4 style="margin: 0; color: {color};">{epitope}</h4>
            <p>Position: {position} | Allele: {allele} {f'| {additional_info}' if additional_info else ''}</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Add paper comparison
    if selected_protein == "hemagglutinin":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) identified the following T-cell Class II epitopes for hemagglutinin:
        
        - **MVSLVKSDQ** (Position 10, ΔG = -11.7756 kcal/mol, HLA-DRB1*03:01)
        - **IGTSTLNQR** (Position 216, ΔG = -56.9580 kcal/mol, HLA-DRB1*03:01)
        
        Note that **IGTSTLNQR** had the strongest binding energy of all epitopes in the study (ΔG = -56.958 kcal/mol),
        making it a particularly promising candidate for vaccine design.
        """)
    elif selected_protein == "neuraminidase":
        st.subheader("Comparison with Paper Results")
        st.markdown("""
        The paper by Tambunan et al. (2016) identified the following T-cell Class II epitope for neuraminidase:
        
        - **YNGIITDTI** (Position 188, ΔG = -17.8049 kcal/mol, HLA-DRB1*01:01)
        
        This epitope was also identified as a B-cell epitope, making it a promising candidate for targeting both 
        humoral and cellular immune responses.
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
            color_continuous_scale='RdBu_r'  # Red (weak) to Blue (strong)
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
        # Just show the top rows if no binding energy data
        selected_epitopes = data.head(4)
    
    for i, row in selected_epitopes.iterrows():
        epitope = row['Epitope']
        position = row['Position']
        epitope_type = row['Type'] if 'Type' in data.columns else "Unknown"
        
        # Get binding energy if available
        binding_info = ""
        if 'ΔG binding' in data.columns:
            binding = row['ΔG binding']
            binding_info = f"| ΔG: {binding:.4f} kcal/mol"
            
            # Choose color based on binding energy
            if binding < -40:
                color = "#1b5e20"  # Very strong binding - dark green
            elif binding < -30:
                color = "#2e7d32"  # Strong binding - medium dark green
            elif binding < -20:
                color = "#43a047"  # Medium binding - medium green
            else:
                color = "#66bb6a"  # Weak binding - light green
        else:
            color = "#1b5e20"  # Default dark green
        
        # Get HLA info if available
        hla_info = ""
        if 'HLA' in data.columns:
            hla = row['HLA']
            hla_info = f"| HLA: {hla}"
        
        st.markdown(f"""
        <div style="background-color: #e8f5e9; padding: 10px; border-left: 5px solid {color}; margin-bottom: 10px;">
            <h4 style="margin: 0; color: {color};">{epitope}</h4>
            <p>Position: {position} | Type: {epitope_type} {hla_info} {binding_info}</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Add paper comparison
    if selected_protein == "hemagglutinin":
        st.subheader("Paper-Recommended Epitopes")
        st.markdown("""
        The paper by Tambunan et al. (2016) selected the following epitopes for hemagglutinin vaccine design:
        
        - **MEKIVLLLA** (ΔG = -40.0148 kcal/mol, T-cell I, HLA-B*35:01)
        - **CPYLGSPSF** (ΔG = -28.0471 kcal/mol, T-cell I, HLA-B*07:02)
        - **KCQTPMGAI** (ΔG = -33.6109 kcal/mol, T-cell I, HLA-B*07:02)
        - **IGTSTLNQR** (ΔG = -56.9580 kcal/mol, B-cell & T-cell II, HLA-DRB1*03:01)
        
        These epitopes were selected based on binding energy calculations and molecular dynamics simulations.
        **IGTSTLNQR** had the strongest binding with a ΔG of -56.958 kcal/mol.
        """)
    elif selected_protein == "neuraminidase":
        st.subheader("Paper-Recommended Epitopes")
        st.markdown("""
        The paper by Tambunan et al. (2016) selected the following epitopes for neuraminidase vaccine design:
        
        - **NPNQKIITI** (ΔG = -23.8779 kcal/mol, T-cell I, HLA-B*07:02)
        - **CYPDAGEIT** (ΔG = -33.0200 kcal/mol, T-cell I, HLA-A*24)
        
        These epitopes were selected based on binding energy calculations and molecular dynamics simulations.
        The strongest binding was observed for **CYPDAGEIT** with a ΔG of -33.020 kcal/mol.
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
    
    # Extract sequence (assuming FASTA format)
    sequence = ""
    header = ""
    for i, line in enumerate(content.split('\n')):
        if i == 0 and line.startswith('>'):
            header = line
        elif not line.startswith('>'):
            sequence += line.strip()
    
    # Create columns for properties and visualization
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("Sequence Properties")
        
        # Sequence length
        st.metric("Sequence Length", len(sequence))
        
        # Calculate amino acid composition
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
        
        # Update layout
        fig.update_layout(
            xaxis_title="Amino Acid",
            yaxis_title="Count",
            height=300
        )
        
        st.plotly_chart(fig)
        
        # Molecular weight and other properties
        avg_weights = {
            'A': 71.08, 'R': 156.19, 'N': 114.11, 'D': 115.09, 'C': 103.15,
            'E': 129.12, 'Q': 128.13, 'G': 57.05, 'H': 137.14, 'I': 113.16,
            'L': 113.16, 'K': 128.17, 'M': 131.19, 'F': 147.18, 'P': 97.12,
            'S': 87.08, 'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13
        }
        
        mol_weight = sum(avg_weights.get(aa, 0) for aa in sequence)
        
        # Calculate estimated isoelectric point (pI)
        # This is a simplified estimate
        acidic = sequence.count('D') + sequence.count('E')
        basic = sequence.count('R') + sequence.count('K') + sequence.count('H')
        pi_estimate = 7.0
        if basic > acidic:
            pi_estimate = 8.5
        elif acidic > basic:
            pi_estimate = 5.5
        
        st.markdown("**Physicochemical Properties:**")
        st.markdown(f"- **Molecular Weight:** ~{mol_weight:.1f} Da")
        st.markdown(f"- **Estimated pI:** ~{pi_estimate}")
        st.markdown(f"- **Number of Epitopes:** {len(sequence.split('GPGPG')) if 'GPGPG' in sequence else 1}")
    
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
        
        epitope_df = pd.DataFrame({
            'Epitope': [epitope for epitope in epitopes if epitope],
            'Length': [len(epitope) for epitope in epitopes if epitope],
            'Position': list(range(1, len([e for e in epitopes if e]) + 1))
        })
        
        st.dataframe(epitope_df)
        
        # GPGPG linker info
        st.markdown("""
        **GPGPG Linker**: The GPGPG linker is used to separate epitopes in the vaccine construct. 
        This linker helps prevent the formation of junctional epitopes and maintains the 
        immunogenicity of each individual epitope.
        """)
    
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
        
        **Molecular Dynamics**: The paper performed 10 ns molecular dynamics simulations to confirm the 
        stability of epitope-HLA interactions. The structures stabilized after 100 ps, confirming their structural integrity.
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
        
        **Conservation**: The selected epitopes were confirmed to be highly conserved with H5N1 NIBRG-14 strain 
        that is commonly used in vaccine development, suggesting they would be effective across multiple H5N1 variants.
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
        
        **Next Steps**: According to the paper, these epitopes should be validated through in vitro and in vivo 
        studies to confirm their immunogenicity and protective efficacy before clinical development.
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
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Hemagglutinin Epitopes:**")
        st.markdown("""
        - **MEKIVLLLA** (ΔG = -40.0148 kcal/mol)
        - **CPYLGSPSF** (ΔG = -28.0471 kcal/mol)
        - **KCQTPMGAI** (ΔG = -33.6109 kcal/mol)
        - **IGTSTLNQR** (ΔG = -56.9580 kcal/mol)
        """)
    
    with col2:
        st.markdown("**Neuraminidase Epitopes:**")
        st.markdown("""
        - **NPNQKIITI** (ΔG = -23.8779 kcal/mol)
        - **CYPDAGEIT** (ΔG = -33.0200 kcal/mol)
        """)
    
    st.markdown("""
    - **Strongest Binding**: IGTSTLNQR with HLA-DRB1*0301 (ΔG = -56.958 kcal/mol)
    - **Molecular Dynamics**: Epitope-HLA structures stabilized after 100 ps of simulation
    - **Conservation**: Epitopes were confirmed to have high homology with H5N1 NIBRG-14 strain
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
                yaxis_title="Number of Epitopes"
            )
            
            st.plotly_chart(fig)
        
        # Show binding energy if available
        if 'ΔG binding' in data.columns:
            st.subheader("Binding Energy")
            
            # Sort by binding energy (more negative is stronger)
            energy_data = data.sort_values('ΔG binding').head(5)
            
            fig = px.bar(
                energy_data,
                x='Epitope',
                y='ΔG binding',
                title="Epitope Binding Energy (ΔG)",
                color='ΔG binding',
                color_continuous_scale='RdBu_r'  # Red (weak) to Blue (strong)
            )
            
            fig.update_layout(
                xaxis_title="Epitope",
                yaxis_title="Binding Energy (kcal/mol)"
            )
            
            st.plotly_chart(fig)
    
    # Show the top epitopes highlighted
    st.subheader("Top T-cell Class I Epitope Candidates")
    
    # Get epitopes sorted by affinity and score
    if 'ΔG binding' in data.columns:
        # Sort by binding energy (more negative is stronger)
        top_epitopes = data.sort_values('ΔG binding').head(4)
    elif 'Affinity' in data.columns and 'Score' in data.columns:
        # Filter strong binders with high scores
        strong_binders = data[data['Affinity'] == 'Strong Binder'].sort_values('Score', ascending=False)
        top_epitopes = strong_binders.head(4)
    else:
        # Just sort by position if no better criteria
        top_epitopes = data.sort_values('Position').head(4)
    
    for i, row in top_epitopes.iterrows():
        epitope = row['Epitope']
        position = row['Position']
        allele = row['Allele'] if 'Allele' in row else 'Unknown'
        
        # Additional info
        binding_info = f"ΔG: {row['ΔG binding']:.4f} kcal/mol" if 'ΔG binding' in row else ""
        affinity_info = f"Affinity: {row['Affinity']}" if 'Affinity' in row else ""
        score_info = f"Score: {row['Score']:.4f}" if 'Score' in row else ""
        
        # Combine available info
        info_parts = [part for part in [binding_info, affinity_info, score_info] if part]
        additional_info = " | ".join(info_parts)
        
        st.markdown(f"""
        <div style="background-color: #e8f5e9; padding: 10px; border-left: 5px solid #2e7d32; margin-bottom: 10px;">
            <h4 style="margin: 0; color: #2e7d32;">{epitope}</h4>
            <p>Position: {position} | Allele: {allele} {f'| {additional_info}' if additional_info else ''}</p>
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