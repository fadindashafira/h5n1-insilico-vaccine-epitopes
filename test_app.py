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
    return Image.open(logo_path)

def run_nextflow_command(cmd, cwd="."):
    """Run a command and capture output in real-time"""
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        text=True,
        cwd=cwd
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
        return None
    
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

def main():
    # Load Kalbe logo
    logo = load_kalbe_logo()
    
    # Show logo in the sidebar
    st.sidebar.image(logo, width=250)
    
    # Main header
    st.markdown('<h1 class="main-header">H5N1 Epitope-Based Vaccine Design</h1>', unsafe_allow_html=True)
    
    # Sidebar for navigation
    st.sidebar.title("Navigation")
    pages = [
        "Home",
        "Configure Pipeline",
        "Run Pipeline",
        "Monitor Progress",
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
            'experiment_id': f'exp_{datetime.now().strftime("%Y%m%d_%H%M%S")}',
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
    
    # Render the selected page
    if selected_page == "Home":
        show_home_page()
    elif selected_page == "Configure Pipeline":
        show_config_page()
    elif selected_page == "Run Pipeline":
        show_run_page()
    elif selected_page == "Monitor Progress":
        show_monitor_page()
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
       - B-cell epitope prediction
       - T-cell Class I epitope prediction
       - T-cell Class II epitope prediction
    3. **Epitope Filtering**: Filter and analyze epitopes based on binding affinity, conservation, and other criteria
    4. **Vaccine Design**: Combine selected epitopes to design optimal vaccine constructs
    5. **Evaluation**: Evaluate the immunogenicity and stability of the designed vaccines
    6. **Optional Molecular Dynamics**: Simulate the stability of epitope-HLA interactions
    
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
        st.metric("Processing Time", "~30 min", "Typical runtime")
    
    # Recent publications
    st.markdown("### Related Research")
    st.markdown("""
    - Tambunan et al. (2016). **Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions**. *Bioinformatics and Biology Insights*, 10, 27-35.
    - Velkov et al. (2013). **The antigenic architecture of the hemagglutinin of influenza H5N1 viruses**. *Molecular Immunology*, 56(4), 705-719.
    - Ahmad et al. (2016). **Immunoinformatic design of a multimeric epitope peptide based vaccine targeting MERS-CoV spike protein: a novel approach**. *Journal of Biomedical Science*, 23(1), 34.
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
        
        bcell_method_options = ["Bepipred-1.0", "Bepipred-2.0", "Chou-Fasman", "Emini", "Karplus-Schulz", "Kolaskar-Tongaonkar", "Parker"]
        bcell_method = st.selectbox(
            "Prediction Method",
            bcell_method_options,
            index=bcell_method_options.index("Bepipred-2.0"),
            help="Algorithm for B-cell epitope prediction"
        )
        
        bcell_threshold = st.slider(
            "Score Threshold",
            0.1, 1.0, 0.5, 0.05,
            help="Minimum score threshold for B-cell epitope prediction"
        )
        
        bcell_length = st.slider(
            "Epitope Length",
            8, 20, 9,
            help="Length of B-cell epitopes to predict"
        )
        
        # T-cell epitope prediction
        st.subheader("T-cell Epitope Prediction")
        
        # Class I HLA
        mhci_alleles = st.text_area(
            "Class I HLA Alleles",
            "HLA-A*02:01, HLA-B*07:02, HLA-B*35:01, HLA-A*11:01",
            help="Comma-separated list of Class I HLA alleles to use for prediction"
        )
        
        mhci_threshold = st.slider(
            "Class I Binding Threshold (nM)",
            50, 1000, 500, 50,
            help="IC50 threshold for MHC Class I binding prediction (lower is more stringent)"
        )
        
        # Class II HLA
        mhcii_alleles = st.text_area(
            "Class II HLA Alleles",
            "HLA-DRB1*03:01, HLA-DRB1*04:01, HLA-DRB1*01:01, HLA-DRB1*07:01",
            help="Comma-separated list of Class II HLA alleles to use for prediction"
        )
        
        mhcii_threshold = st.slider(
            "Class II Binding Threshold (nM)",
            50, 1000, 500, 50,
            help="IC50 threshold for MHC Class II binding prediction (lower is more stringent)"
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
            5, 30, 10,
            help="Maximum number of epitopes to include in the vaccine construct"
        )
        
        # These advanced options would be passed to the pipeline
        # But for simplicity we're not implementing all param passing now
        st.info("These advanced options can be configured for specialized applications. Default values are suitable for most cases.")

def show_run_page():
    st.markdown('<h2 class="sub-header">Run Pipeline</h2>', unsafe_allow_html=True)
    
    # Check if a pipeline is already running
    if st.session_state.run_status.get('is_running', False):
        st.markdown("""
        <div class="warning-box">
        A pipeline execution is already in progress. Please wait for it to complete or check the 
        Monitor Progress page for updates.
        </div>
        """, unsafe_allow_html=True)
        
        # Show current run info
        st.subheader("Current Run Information")
        st.write(f"Run ID: {st.session_state.run_status['run_id']}")
        st.write(f"Started at: {st.session_state.run_status['start_time']}")
        
        # Option to stop the pipeline
        if st.button("Stop Pipeline"):
            if st.session_state.run_status['process']:
                try:
                    # Terminate the process
                    st.session_state.run_status['process'].terminate()
                    st.session_state.run_status['is_running'] = False
                    st.success("Pipeline stopped successfully.")
                except Exception as e:
                    st.error(f"Error stopping pipeline: {str(e)}")
            else:
                st.session_state.run_status['is_running'] = False
                st.success("Pipeline status reset.")
        
        st.markdown("Go to the **Monitor Progress** page to see detailed status updates.")
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
    
    # Run button
    if st.button("Run Pipeline", key="run_pipeline_button"):
        # Create a unique run ID
        run_id = f"run_{int(time.time())}"
        
        # Update run status
        st.session_state.run_status['is_running'] = True
        st.session_state.run_status['run_id'] = run_id
        st.session_state.run_status['start_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        st.session_state.run_status['experiment_output'] = experiment_output
        
        try:
            # Execute the command
            process = run_nextflow_command(cmd)
            st.session_state.run_status['process'] = process
            
            # Set the log file location
            log_file = ".nextflow.log"
            st.session_state.run_status['log_file'] = log_file
            
            st.success(f"Pipeline started successfully! Run ID: {run_id}")
            st.info("Go to the **Monitor Progress** page to see status updates.")
        except Exception as e:
            st.error(f"Error starting pipeline: {str(e)}")
            st.session_state.run_status['is_running'] = False

def show_monitor_page():
    st.markdown('<h2 class="sub-header">Monitor Pipeline Progress</h2>', unsafe_allow_html=True)
    
    # Check if a pipeline is running
    if not st.session_state.run_status.get('is_running', False):
        st.markdown("""
        <div class="info-box">
        No pipeline is currently running. Go to the Run Pipeline page to start a new execution.
        </div>
        """, unsafe_allow_html=True)
        
        # Allow checking previous runs
        st.markdown("### Check Previous Runs")
        
        # Find all nextflow log files
        log_files = sorted(glob.glob(".nextflow.log*"))
        
        if log_files:
            selected_log = st.selectbox("Select a log file to view", log_files)
            
            if selected_log:
                # Parse and display the log
                log_info = parse_nextflow_log(selected_log)
                
                if log_info:
                    st.markdown(f"**Status:** {log_info.get('status', 'Unknown')}")
                    
                    if 'session_id' in log_info:
                        st.markdown(f"**Session ID:** {log_info['session_id']}")
                    
                    if 'command' in log_info:
                        st.markdown("**Command:**")
                        st.code(log_info['command'])
                    
                    if 'processes' in log_info and log_info['processes']:
                        st.markdown("**Process Status:**")
                        process_chart = create_process_status_chart(log_info['processes'])
                        if process_chart:
                            st.plotly_chart(process_chart)
                
                # Display the raw log
                with st.expander("View Raw Log"):
                    try:
                        with open(selected_log, 'r') as f:
                            log_content = f.read()
                        st.markdown(f'<div class="terminal">{log_content}</div>', unsafe_allow_html=True)
                    except Exception as e:
                        st.error(f"Error reading log file: {str(e)}")
        else:
            st.info("No Nextflow log files found.")
        
        return
    
    # Display current run information
    st.markdown("### Current Run Information")
    st.write(f"Run ID: {st.session_state.run_status['run_id']}")
    st.write(f"Started at: {st.session_state.run_status['start_time']}")
    st.write(f"Experiment output: {st.session_state.run_status['experiment_output']}")
    
    # Auto-refresh option
    auto_refresh = st.checkbox("Auto-refresh (every 10 seconds)", value=True)
    
    # Create columns for status and controls
    col1, col2 = st.columns([3, 1])
    
    with col2:
        # Stop button
        if st.button("Stop Pipeline"):
            if st.session_state.run_status['process']:
                try:
                    # Terminate the process
                    st.session_state.run_status['process'].terminate()
                    st.session_state.run_status['is_running'] = False
                    st.success("Pipeline stopped successfully.")
                except Exception as e:
                    st.error(f"Error stopping pipeline: {str(e)}")
            else:
                st.session_state.run_status['is_running'] = False
                st.success("Pipeline status reset.")
    
    # Check process status and update if needed
    process = st.session_state.run_status.get('process')
    if process:
        # Check if process has completed
        return_code = process.poll()
        
        if return_code is not None:
            # Process has completed
            st.session_state.run_status['is_running'] = False
            
            if return_code == 0:
                st.success("Pipeline completed successfully!")
            else:
                st.error(f"Pipeline exited with error code: {return_code}")
                
                # Display error output
                stderr = process.stderr.read()
                if stderr:
                    st.markdown("**Error output:**")
                    st.markdown(f'<div class="terminal">{stderr}</div>', unsafe_allow_html=True)
    
    # Read and parse the log file
    log_file = st.session_state.run_status.get('log_file', '.nextflow.log')
    log_info = parse_nextflow_log(log_file)
    
    if log_info:
        with col1:
            # Display status
            status = log_info.get('status', 'RUNNING')
            if status == 'COMPLETED':
                st.success(f"Status: {status}")
            elif status == 'ERROR':
                st.error(f"Status: {status}")
            elif status == 'CANCELLED':
                st.warning(f"Status: {status}")
            else:
                st.info(f"Status: {status}")
        
        # Create process status visualization
        if 'processes' in log_info and log_info['processes']:
            st.subheader("Process Status")
            process_chart = create_process_status_chart(log_info['processes'])
            if process_chart:
                st.plotly_chart(process_chart)
    
    # Find and display timeline if available
    experiment_output = st.session_state.run_status.get('experiment_output')
    if experiment_output:
        timeline_file = f"{experiment_output}/reports/timeline.html"
        if os.path.exists(timeline_file):
            st.subheader("Pipeline Timeline")
            st.markdown(f"Timeline report available at: `{timeline_file}`")
            
            # Parse timeline data
            timeline_data = parse_timeline_file(timeline_file)
            if timeline_data:
                # Display summary of timeline data
                st.markdown("**Timeline Summary:**")
                
                # Calculate some stats
                tasks = len(timeline_data)
                completed = sum(1 for task in timeline_data if task.get('status') == 'COMPLETED')
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total Tasks", tasks)
                with col2:
                    st.metric("Completed", completed)
                with col3:
                    st.metric("Completion", f"{int(completed/tasks*100)}%" if tasks > 0 else "0%")
                
                # Offer to open the timeline in a new tab
                st.markdown(f'<a href="{timeline_file}" target="_blank">Open Timeline Visualization</a>', unsafe_allow_html=True)
    
    # Display log content
    st.subheader("Log Output")
    log_container = st.empty()
    
    # Read and display the log file
    try:
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                log_content = f.read()
            log_container.markdown(f'<div class="terminal">{log_content}</div>', unsafe_allow_html=True)
        else:
            log_container.warning(f"Log file not found: {log_file}")
    except Exception as e:
        log_container.error(f"Error reading log file: {str(e)}")
    
    # Auto-refresh logic
    if auto_refresh and st.session_state.run_status.get('is_running', False):
        st.markdown("Page will refresh automatically every 10 seconds...")
        st.markdown(f'<meta http-equiv="refresh" content="10">', unsafe_allow_html=True)

def show_results_page():
    st.markdown('<h2 class="sub-header">View Results</h2>', unsafe_allow_html=True)
    
    # Find all experiment directories
    result_dirs = glob.glob("results/results_*")
    
    if not result_dirs:
        st.markdown("""
        <div class="info-box">
        No results found. Please run the pipeline first to generate results.
        </div>
        """, unsafe_allow_html=True)
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
        return
    
    # Create tabs for different result types
    tabs = ["B-cell Epitopes", "T-cell Class I Epitopes", "T-cell Class II Epitopes", "Combined Epitopes", "Vaccine Constructs"]
    selected_tab = st.radio("Result Type", tabs)
    
    # Display the selected result type
    if selected_tab == "B-cell Epitopes":
        show_bcell_results(results)
    elif selected_tab == "T-cell Class I Epitopes":
        show_tcell_i_results(results)
    elif selected_tab == "T-cell Class II Epitopes":
        show_tcell_ii_results(results)
    elif selected_tab == "Combined Epitopes":
        show_combined_results(results)
    elif selected_tab == "Vaccine Constructs":
        show_vaccine_results(results)
    
    # Show reports if available
    st.markdown("### Pipeline Reports")
    reports_dir = f"{selected_experiment}/reports"
    
    if os.path.exists(reports_dir):
        report_files = glob.glob(f"{reports_dir}/*")
        
        if report_files:
            for report in report_files:
                report_name = os.path.basename(report)
                st.markdown(f"- [{report_name}]({report})")
        else:
            st.info(f"No report files found in {reports_dir}")
    else:
        st.info(f"Reports directory not found: {reports_dir}")

def show_bcell_results(results):
    st.markdown("### B-cell Epitope Prediction Results")
    
    if not results['bcell_epitopes']:
        st.info("No B-cell epitope results found.")
        return
    
    # Create a dropdown to select protein
    proteins = [item['protein'] for item in results['bcell_epitopes']]
    selected_protein = st.selectbox("Select Protein", proteins)
    
    # Find the selected protein's data
    for item in results['bcell_epitopes']:
        if item['protein'] == selected_protein:
            data = item['data']
            
            # Display the data table
            st.subheader(f"{selected_protein} B-cell Epitopes")
            st.dataframe(data)
            
            # Create a bar chart of epitope scores
            if 'Score' in data.columns:
                st.subheader("Epitope Scores")
                fig = px.bar(
                    data,
                    x='Epitope' if 'Epitope' in data.columns else data.index,
                    y='Score',
                    title=f"{selected_protein} B-cell Epitope Scores",
                    color='Score',
                    color_continuous_scale='greens'
                )
                st.plotly_chart(fig)
            
            # Provide download link
            st.download_button(
                label=f"Download {selected_protein} B-cell epitopes",
                data=data.to_csv(index=False),
                file_name=f"{selected_protein}_bcell_epitopes.csv",
                mime="text/csv"
            )
            
            break

def show_tcell_i_results(results):
    st.markdown("### T-cell Class I Epitope Prediction Results")
    
    if not results['tcell_i_epitopes']:
        st.info("No T-cell Class I epitope results found.")
        return
    
    # Create a dropdown to select protein
    proteins = [item['protein'] for item in results['tcell_i_epitopes']]
    selected_protein = st.selectbox("Select Protein", proteins)
    
    # Find the selected protein's data
    for item in results['tcell_i_epitopes']:
        if item['protein'] == selected_protein:
            data = item['data']
            
            # Display the data table
            st.subheader(f"{selected_protein} T-cell Class I Epitopes")
            st.dataframe(data)
            
            # Create visualizations based on available columns
            if 'Allele' in data.columns and 'Score' in data.columns:
                st.subheader("Epitope Binding by Allele")
                
                # Group by allele
                allele_data = data.groupby('Allele')['Score'].mean().reset_index()
                
                fig = px.bar(
                    allele_data,
                    x='Allele',
                    y='Score',
                    title=f"{selected_protein} T-cell Class I Epitope Binding by Allele",
                    color='Score',
                    color_continuous_scale='greens'
                )
                st.plotly_chart(fig)
            
            # Provide download link
            st.download_button(
                label=f"Download {selected_protein} T-cell Class I epitopes",
                data=data.to_csv(index=False),
                file_name=f"{selected_protein}_tcell_i_epitopes.csv",
                mime="text/csv"
            )
            
            break

def show_tcell_ii_results(results):
    st.markdown("### T-cell Class II Epitope Prediction Results")
    
    if not results['tcell_ii_epitopes']:
        st.info("No T-cell Class II epitope results found.")
        return
    
    # Create a dropdown to select protein
    proteins = [item['protein'] for item in results['tcell_ii_epitopes']]
    selected_protein = st.selectbox("Select Protein", proteins)
    
    # Find the selected protein's data
    for item in results['tcell_ii_epitopes']:
        if item['protein'] == selected_protein:
            data = item['data']
            
            # Display the data table
            st.subheader(f"{selected_protein} T-cell Class II Epitopes")
            st.dataframe(data)
            
            # Create visualizations based on available columns
            if 'Allele' in data.columns and 'Score' in data.columns:
                st.subheader("Epitope Binding by Allele")
                
                # Group by allele
                allele_data = data.groupby('Allele')['Score'].mean().reset_index()
                
                fig = px.bar(
                    allele_data,
                    x='Allele',
                    y='Score',
                    title=f"{selected_protein} T-cell Class II Epitope Binding by Allele",
                    color='Score',
                    color_continuous_scale='greens'
                )
                st.plotly_chart(fig)
            
            # Provide download link
            st.download_button(
                label=f"Download {selected_protein} T-cell Class II epitopes",
                data=data.to_csv(index=False),
                file_name=f"{selected_protein}_tcell_ii_epitopes.csv",
                mime="text/csv"
            )
            
            break

def show_combined_results(results):
    st.markdown("### Combined Epitopes Results")
    
    if not results['combined_epitopes']:
        st.info("No combined epitope results found.")
        return
    
    # Create a dropdown to select protein
    proteins = [item['protein'] for item in results['combined_epitopes']]
    selected_protein = st.selectbox("Select Protein", proteins)
    
    # Find the selected protein's data
    for item in results['combined_epitopes']:
        if item['protein'] == selected_protein:
            data = item['data']
            
            # Display the data table
            st.subheader(f"{selected_protein} Combined Epitopes")
            st.dataframe(data)
            
            # Create visualizations based on available columns
            if 'Type' in data.columns and 'Score' in data.columns:
                st.subheader("Epitopes by Type")
                
                # Group by type
                type_counts = data['Type'].value_counts().reset_index()
                type_counts.columns = ['Type', 'Count']
                
                fig = px.pie(
                    type_counts,
                    names='Type',
                    values='Count',
                    title=f"{selected_protein} Epitope Types",
                    color_discrete_sequence=px.colors.sequential.Greens
                )
                st.plotly_chart(fig)
            
            # Provide download link
            st.download_button(
                label=f"Download {selected_protein} combined epitopes",
                data=data.to_csv(index=False),
                file_name=f"{selected_protein}_combined_epitopes.csv",
                mime="text/csv"
            )
            
            break

def show_vaccine_results(results):
    st.markdown("### Vaccine Constructs")
    
    if not results['vaccine_constructs']:
        st.info("No vaccine construct results found.")
        return
    
    # Create a dropdown to select protein
    proteins = [item['protein'] for item in results['vaccine_constructs']]
    selected_protein = st.selectbox("Select Protein", proteins)
    
    # Find the selected protein's data
    for item in results['vaccine_constructs']:
        if item['protein'] == selected_protein:
            content = item['content']
            
            # Display the FASTA content
            st.subheader(f"{selected_protein} Vaccine Construct")
            st.text(content)
            
            # Extract the sequence (assuming FASTA format)
            sequence = ""
            for line in content.split('\n'):
                if not line.startswith('>'):
                    sequence += line.strip()
            
            # Display sequence properties
            st.subheader("Sequence Properties")
            
            length = len(sequence)
            
            # Calculate amino acid composition
            aa_count = {}
            for aa in sorted(set(sequence)):
                aa_count[aa] = sequence.count(aa)
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.metric("Sequence Length", length)
                
                # Display amino acid composition
                st.markdown("**Amino Acid Composition:**")
                
                aa_df = pd.DataFrame({
                    'Amino Acid': list(aa_count.keys()),
                    'Count': list(aa_count.values())
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
                # Create a visual representation of the vaccine construct
                # This is a simplified visualization - in a real app you'd want more sophisticated visualization
                st.markdown("**Vaccine Construct Visualization:**")
                
                # Detect linkers (assuming GPGPG)
                parts = sequence.split("GPGPG")
                
                # Create a simple visualization with colored segments for epitopes and linkers
                fig = go.Figure()
                
                y_pos = 1
                x_start = 0
                
                for i, part in enumerate(parts):
                    # Add epitope segment
                    if part:
                        fig.add_trace(go.Scatter(
                            x=[x_start, x_start + len(part)],
                            y=[y_pos, y_pos],
                            mode='lines',
                            line=dict(width=10, color='#43a047'),
                            name=f'Epitope {i+1}'
                        ))
                        
                        # Add text label
                        fig.add_annotation(
                            x=(x_start + x_start + len(part)) / 2,
                            y=y_pos + 0.1,
                            text=f'E{i+1}',
                            showarrow=False
                        )
                        
                        x_start += len(part)
                    
                    # Add linker if not the last segment
                    if i < len(parts) - 1:
                        fig.add_trace(go.Scatter(
                            x=[x_start, x_start + 5],  # 5 for "GPGPG"
                            y=[y_pos, y_pos],
                            mode='lines',
                            line=dict(width=10, color='#ffb300'),
                            name='Linker'
                        ))
                        
                        # Add text label for linker
                        fig.add_annotation(
                            x=x_start + 2.5,
                            y=y_pos + 0.1,
                            text='L',
                            showarrow=False
                        )
                        
                        x_start += 5
                
                # Update layout
                fig.update_layout(
                    title='Vaccine Construct Structure',
                    xaxis=dict(title='Position'),
                    showlegend=False,
                    height=200,
                    margin=dict(l=20, r=20, t=40, b=20)
                )
                
                st.plotly_chart(fig)
            
            # Provide download link
            st.download_button(
                label=f"Download {selected_protein} vaccine construct",
                data=content,
                file_name=f"{selected_protein}_vaccine_construct.fasta",
                mime="text/plain"
            )
            
            break

def show_about_page():
    st.markdown('<h2 class="sub-header">About</h2>', unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
    <h3>H5N1 Epitope-Based Vaccine Design Pipeline</h3>
    <p>This application provides a user-friendly interface for designing epitope-based vaccines against H5N1 avian influenza virus using immunoinformatics approaches.</p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    ### Background
    
    Avian influenza, caused by H5N1 virus, represents a significant public health concern due to its high mortality rate and pandemic potential. Developing effective vaccines against H5N1 is crucial for pandemic preparedness.
    
    This tool implements the methodology described in the paper by Tambunan et al. (2016), "Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions," to predict potential epitopes from the hemagglutinin and neuraminidase proteins of H5N1 virus.
    
    ### Pipeline Overview
    
    The pipeline follows these key steps:
    
    1. **Sequence Retrieval**: Obtains hemagglutinin and neuraminidase protein sequences from NCBI.
    2. **Epitope Prediction**:
       - Predicts B-cell epitopes using algorithms like Bepipred
       - Predicts T-cell epitopes for Class I and Class II HLA alleles
    3. **Epitope Filtering**: Selects high-quality epitopes based on binding affinity, conservation, and other criteria.
    4. **Vaccine Design**: Combines selected epitopes with appropriate linkers to create vaccine constructs.
    5. **Evaluation**: Assesses the immunogenicity and stability of the designed vaccines.
    
    ### Implementation
    
    The pipeline is implemented using Nextflow, a workflow management system that enables scalable and reproducible scientific workflows. The application integrates with a suite of bioinformatics tools for epitope prediction, structural modeling, and analysis.
    
    ### Credits
    
    **Pipeline Development**: Kalbe Bioinformatics Team
    
    **Core Methodology**:
    - Tambunan et al. (2016). Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions. *Bioinformatics and Biology Insights*, 10, 27-35.
    
    **User Interface**: Streamlit-based application
    
    ### License
    
    This software is provided for research and educational purposes only. Not for clinical use.
    """)
    
    # Version information
    st.markdown("### Version Information")
    st.markdown("- **Application Version**: 1.0.0")
    st.markdown("- **Pipeline Version**: 1.0.0")
    st.markdown("- **Last Updated**: April 2025")

# Run the main function
if __name__ == "__main__":
    main()