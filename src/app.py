import streamlit as st
import os
import pandas as pd
import matplotlib.pyplot as plt

from scan_p53 import reader, get_allPWM, calculate_background, read_input, to_CSV
from plots_hit import plot_motif_hits

st.set_page_config(page_title="p53 Motif Scanner", layout="wide")
st.title("p53 Motif Scanner")
st.markdown("Upload a DNA sequence to scan for p53 binding sites using a statistical Position Weight Matrix (PWM).")

with st.sidebar:
    st.header("Settings")
    tf_name = st.text_input("Transcription Factor", value="TP53")
    p_value = st.number_input("P-value Cutoff (Significance)", 
                              value=0.0001, 
                              min_value=1e-10, 
                              max_value=0.05, 
                              format="%.6f",
                              help="Lower value = Stricter/Stronger sites only")
    input_method = st.radio("Input Method", ["Paste Text", "Upload FASTA"])

seq_list = []

if input_method == "Paste Text":
    raw_text = st.text_area("Paste DNA Sequence here (FASTA format or raw seq):", height=200)
    if raw_text:
        if ">" in raw_text:
            with open("temp_input.fasta", "w") as f:
                f.write(raw_text)
            seq_list = read_input("temp_input.fasta", is_file=True)
        else:
            seq_list = read_input(raw_text, is_file=False)

elif input_method == "Upload FASTA":
    uploaded_file = st.file_uploader("Choose a FASTA file", type=["fasta", "txt"])
    if uploaded_file is not None:
        string_data = uploaded_file.getvalue().decode("utf-8")
        with open("temp_input.fasta", "w") as f:
            f.write(string_data)
        seq_list = read_input("temp_input.fasta", is_file=True)

if st.button("Scan Sequence"):
    if not seq_list:
        st.error("Please provide DNA sequence data first.")
    else:
        with st.spinner(f'Fetching Motifs {tf_name} with p_value = {p_value} and Scanning...'):
            try:
                ids = reader(tf_name)
                gc = calculate_background(seq_list)
                pwm_all = get_allPWM(ids, gc)
                
                output_csv = "temp_results.csv"
                to_CSV(seq_list, pwm_all, ids, gc, filename=output_csv, pval_cutoff=p_value)
                
                st.success("Scan Complete!")
                if os.path.exists(output_csv) and os.path.getsize(output_csv) > 0:
                    df = pd.read_csv(output_csv)
                   
                    
                    if not df.empty:
                        st.session_state['scan_results'] = df
                        st.session_state['unique_seqs'] = df["Sequence_ID"].unique()
                        st.success("Scan Complete! Results stored.")
                    else:
                        st.warning(f"No binding sites found for {tf_name} with P < {p_value}")
                        st.session_state['scan_results'] = pd.DataFrame() 
                else:
                    st.error("No results generated.")
                    
            except Exception as e:
                st.error(f"An error occurred: {e}")
                
scan_results = st.session_state.get('scan_results')
if st.session_state['scan_results'] is not None and not st.session_state['scan_results'].empty:
    
    df = st.session_state['scan_results']
    unique_seqs = st.session_state['unique_seqs']
    n_seqs = len(unique_seqs)

    st.divider()
    st.success(f"Loaded {len(df)} binding sites across {n_seqs} sequences.")
    
    st.subheader("Visualization")
    selected_seq = st.selectbox("Choose a Sequence to Inspect:", unique_seqs)
    
    if selected_seq:
        subset_df = df[df["Sequence_ID"] == selected_seq]
        try:
            fig = plot_motif_hits(subset_df, seq_name=selected_seq)
            st.pyplot(fig)
        except Exception as e:
            st.error(f"Error plotting sequence: {e}")
        
        st.divider()
        col1, col2 = st.columns(2)
        with col1:
            st.write(f"**Hits for {selected_seq}**")
            st.dataframe(subset_df)
        
        with col2:
            st.write("**All Hits (Summary)**")
            st.dataframe(df)
    csv_data = df.to_csv(index=False).encode('utf-8')
    st.download_button("Download All Results (CSV)", csv_data, f"{tf_name}_hits.csv", "text/csv")

