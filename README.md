# Universal DNA Motif Scanner

Project Overview: A statistical bioinformatics web application that identifies and visualizes transcription factor binding sites in DNA sequences.

Core Logic: Uses Position Weight Matrices (PWMs) from the JASPAR database to calculate binding affinity.

Statistical Rigor: Implements Monte Carlo simulations to calculate P-values and applies dynamic GC-content background correction to eliminate false positives.

Tech Stack: Python, Streamlit, Biopython, Matplotlib, Pandas.

 
Case Study: The p53 Tumor Suppressor
The scanner was tasked with finding p53 binding sites in the human CDKN1A (p21) promoter region (GenBank U24170.1).

<img width="1460" alt="p53 Binding Landscape" src="https://github.com/user-attachments/assets/b827a40e-f66c-46cc-b813-8ceae793aee2" />

The Signal (Red Spike): The algorithm correctly identified the canonical p53 Response Element (RE-1) at the clinically verified location (~160bp upstream) with a high Log-Odds score (>2.0).

The Noise (Grey Bars): The "Half-Site" markers indicate weak, transient binding events. The clear separation between the Red Signal and Grey Noise demonstrates the tool's high Signal-to-Noise Ratio (SNR).


