import matplotlib.pyplot as plt
import csv
import pandas as pd
import argparse

def plot_motif_hits(df, seq_name = "Sequence"):
    full_pos = []
    full_score = []
    half_pos  = []
    half_score = []

    for ind, row in df.iterows():
        pos = int(row["Position"])
        score = float(row["Score"])
        site = int(row["Full/ Half-site"])
        if (score > 0):
            if site == 1:
                full_pos.append(pos)
                full_score.append(score)
            else:
                half_pos.append(pos)
                half_score.append(score)
    fig, ax =  plt.subplots(figsize=(12,8))
    ax.axhline(0, color = "black", linewidth = 1)
    if full_pos:
        ax.bar(full_pos, 
                full_score, 
                color = "red", 
                width = 10, 
                label = "Full Site (Strong)",
                zorder = 10)
    if half_pos:
        ax.bar(half_pos, half_score, 
                color = "grey",
                width = 10, 
                alpha = 0.3, 
                label = "Half Site (Weak)",
                zorder = 1)
    ax.set_title(f"Binding Sites on {seq_name}")
    ax.set_xlabel("Position in Sequence (bp)")
    ax.set_ylabel("Binding Score (Log-Odds Score)")
    ax.legend(loc = 'upper right')
    ax.grid(axis= 'y', linestyle = "--", alpha = 0.5)
    return fig

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", default="../output/p21_results.csv", help="Path to CSV results")
    args = parser.parse_args()
    
    plot_motif_hits(args.file)

    plt.show()
