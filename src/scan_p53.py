from Bio import motifs
import requests
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import motifs
import os
import math
import csv
import random
import bisect
import argparse
import sys

bases = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
background = 0.25

def reader(tf_name = "TP53"):
    url = f"https://jaspar.elixir.no/api/v1/matrix/?search={tf_name}&collection=CORE&tax_group=vertebrates"
    response = requests.get(url)
    res = response.json()["results"]
    ids = []
    for motif in res:
        props = (motif["matrix_id"], motif["url"])
        ids.append(props)   
    return ids

def read_fasta(file_path):
    sequences = []  
    for line in SeqIO.parse(file_path, "fasta"):
        sequences.append((line.id, str(line.seq).upper()))
    return sequences

def read_input(input_data, is_file=False):
    sequences = []

    if is_file:
        return read_fasta(input_data)
    else:
        for idx, line in enumerate(input_data.splitlines(), start=0):
            line = line.strip().upper()
            if len(line) > 0:
                sequences.append((f"seq{idx}", line))

    return sequences

def get_PWM(pfm, pseudocount = 0.1):
    pwm = []
    length = len(pfm["A"])

    for i in range(length):
        column = []
        total = sum(pfm[b][i] for b, pos in bases.items()) + pseudocount*4
        for b in bases:
            freq = (pfm[b][i] + pseudocount)/ total
            score = math.log2(freq/ background) 
            column.append(score)
        pwm.append(column)
    return pwm

def get_allPWM(ids):
    pwm_all = []
    for props in ids:
        id = props[0]
        url = props[1]
        response = requests.get(url)
        pwm = response.json()["pfm"]
        pwm_all.append(get_PWM(pwm))
    return pwm_all

def scan(seq, pwm, sign, threshold_fraction = 0.5, fixed_threshold= None):
    hits = []
    seq = seq.upper()
    motif_length = len(pwm)
    seq_length = len(seq)
    max_per_base = get_maxscore_perbase(pwm)
    threshold = fixed_threshold if fixed_threshold is not None else max_per_base * threshold_fraction

    for i in range(seq_length - motif_length + 1):
        motif = seq[i:i + motif_length]
        score = 0
        flag = True
        for j, base in enumerate(motif):
            if base in bases:
                score += pwm[j][bases[base]]
            else:
                flag = False
                break
        if flag:
            avg_score = score/motif_length
            if avg_score >= threshold:
                    if (sign == "+"):
                        hits.append((i, sign, round(avg_score,3)))
                    else:
                        rev_pos = seq_length - motif_length - i
                        hits.append((rev_pos, sign,round(avg_score,3)))

    return hits

def get_maxscore(pwm):
    max_score = 0
    for col in pwm:
        max_score += max(col)
    return max_score

def get_maxscore_perbase(pwm):
    max_score = 0
    for col in pwm:
        max_score += max(col)
    return max_score/len(pwm)

def reverse_complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(complement.get(b,b) for b in seq[::-1])

def generate_null_dist(pwm, n = 10000):
    n_scores = []
    motif_length = len(pwm)
    options = [0,1,2,3]

    for _ in range(n):
        random_seq = [random.choice(options) for _ in range(motif_length)]
        score = 0
        for pos, base_ind in enumerate(random_seq):
            score += pwm[pos][base_ind]
        n_scores.append(score)
    n_scores.sort()
    return n_scores

def p_val(score, n_score):
    n = len(n_score)
    ind = bisect.bisect_left(n_score, score)
    num_better = n - ind
    return (num_better + 1)/ (n + 1) #pseudocount = 1 

def threshold_p(p, n_scores):
    n = len(n_scores)
    num_top = int(p*n)
    if num_top < 1:
        num_top = 1
    ind = n - num_top
    print(n_scores[ind])
    return n_scores[ind]

def to_CSV(seq_list, pwm_all, motif_ids, filename, threshold_fraction= 0.5, pval_cutoff=0.001):
    directory = os.path.dirname(filename)
    if directory: 
        os.makedirs(directory, exist_ok=True)
    with open(filename, "w", newline="") as f:
        write = csv.writer(f)
        write.writerow(["Sequence_ID", "Motif_ID", "Position", "Direction", "Score", "P-value", "Full/ Half-site"])
        for idx, pwm in enumerate(pwm_all):
            hs_len = len(pwm)//2
            null_dist_full = generate_null_dist(pwm)
            null_dist_half = generate_null_dist(pwm[:hs_len])
            th_full = threshold_p(pval_cutoff, null_dist_full)
            th_half = threshold_p(pval_cutoff, null_dist_half)
            for seq_id, seq in seq_list:
                full_site = pwm
                half_site = pwm[:hs_len]

                hits_full = scan(seq, full_site, '+', threshold_fraction)
                hits_full += scan(reverse_complement(seq), full_site, '-', fixed_threshold= th_full)

                hits_half = scan(seq, half_site, '+', threshold_fraction)
                hits_half += scan(reverse_complement(seq), half_site, '-', fixed_threshold= th_half)
                
                hit_groups = [
                    (hits_full, 1, null_dist_full), 
                    (hits_half, 0, null_dist_half)
                ]

                for hits, site_type, n_dist in hit_groups:
                    hits_sorted = sorted(hits, key=lambda x: x[2], reverse=True)
                    for pos, sign, score in hits_sorted:
                        p_value = p_val(score, n_dist)
                        write.writerow([seq_id, motif_ids[idx][0], pos, sign, score, p_value, site_type])

def main():
    parser = argparse.ArgumentParser(description="Scan DNA sequences for TF binding sites using statistical thresholding.")
    parser.add_argument("-m", "--motif", type=str, required=True, help="Transcription Factor name (e.g., TP53, MA0106.1)")
    parser.add_argument("-f", "--fasta", type=str, default="../tests/test_seqs.fasta", help="Path to input FASTA")
    parser.add_argument("-o", "--out", type=str, default="../output/motif_hits.csv", help="Output filename")
    parser.add_argument("-p", "--pvalue", type=float, default=0.001, help="P-value cutoff (default: 0.001)")
    args = parser.parse_args()
    print(f"Scan for {args.motif}")
    print(f"Fetching matrix for {args.motif} from JASPAR...")
    motif_ids = reader(args.motif)
    
    if not motif_ids:
        print(f"Error: No motifs found for '{args.motif}'")
        sys.exit(1)
        
    pwm_all = get_allPWM(motif_ids)

    print(f"Reading sequences from {args.fasta}...")
    if not os.path.exists(args.fasta):
        print(f"Error: File '{args.fasta}' not found.")
        sys.exit(1)
        
    seq_list = read_input(args.fasta, is_file=True)
    pwm_all = get_allPWM(motif_ids)
    to_CSV(seq_list, pwm_all, motif_ids, filename=args.out, pval_cutoff=args.pvalue)

if __name__ == "__main__":
    main()

