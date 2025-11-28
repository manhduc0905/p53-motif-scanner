from Bio import Entrez, SeqIO 
Entrez.email = "bioinformatics@edu.org"

def fetch_p21_promoter():
    handle = Entrez.efetch(db = "nucleotide", id="NG_016708.1", rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    promoter_seq = record.seq
    output_path = "../tests/p53_seq.fasta"
    with open(output_path, "w") as f:
        f.write(">p53\n")
        f.write(str(promoter_seq))
if __name__ == "__main__":
    fetch_p21_promoter()
