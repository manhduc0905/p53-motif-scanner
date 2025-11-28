import streamlit as st
import os
import pandas as pd
import matplotlib.pyplot as plt

from scan_p53 import reader, get_allPWM, calculate_background, read_input, to_CSV
from plots_hit import plot_motif_hits

st.set_page_config(page_title="Motif Scanner", layout="wide")
st.title("Motif Scanner")
st.markdown("Upload a DNA sequence to scan for a choosen transcription factor binding sites using Position Weight Matrix (PWM).")

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
                    st.dataframe(df)
                    
                    if not df.empty:
                        st.subheader("Visualization")
                        fig = plot_motif_hits(output_csv)
                        st.pyplot(fig)
                    else:
                        st.warning("No hits found above the threshold.")
                    with open(output_csv, "rb") as file:
                        st.download_button("Download CSV", file, "motif_results.csv")
                else:
                    st.error("No results generated.")
                    
            except Exception as e:
                st.error(f"An error occurred: {e}")
