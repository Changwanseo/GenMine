from lib import GenMiner as Gen
import argparse
from datetime import datetime, date

import os

# Default arguments
# email = ""
# genus_term = "" #Genus

date_start = -1
date_end = -1

# Argument input for running in commandline mode
parser = argparse.ArgumentParser(description="Download genbank records")
parser.add_argument("--email", "-e", help="entrez email", type=str)
parser.add_argument("--genus", "-g", nargs="*", help="genus to find", type=str)
parser.add_argument("--additional", "-a", nargs="*", help="additional terms", type=str)
parser.add_argument("--accession", "-c", help="Accession list file to get", type=str)
parser.add_argument("--out", "-o", help="out path", type=str)
parser.add_argument(
    "--max",
    "-m",
    help="maximum length of sequence to search, in order to remove genomic sequences",
    type=int,
)
parser.add_argument(
    "--start",
    "-s",
    help="starting date of finding records, as term of YYYY-MM-DD",
    type=str,
)
parser.add_argument(
    "--end", "-f", help="final date of finding records, as term of YYYY-MM-DD", type=str
)
# parser.add_argument("--start", "-d1", help = "starting search date", type=int)
# parser.add_argument("--end", "-d2", help = "starting end date", type=int)

args = parser.parse_args()

if args.email != None:
    email = args.email
else:
    print("No emails accepted, this might be bad for NCBI connection")

genus_term = None
accession_file = None

if args.accession != None:
    accession_file = args.accession
elif args.genus != None:
    genus_term = args.genus
else:
    print("No genus term nor accession numbers inserted. Aborted")
    raise Exception

if args.additional != None:
    additional_term = list(args.additional)

if args.out != None:
    path_out = args.out

if args.max != None:
    max_len = args.max

if date_start != None:
    date_start = args.start

path_localgb = f"{path_out}.json"
path_localgb_xlsx = f"{path_out}.xlsx"

date = date.today().strftime("%Y-%m-%d")
path_out = f"{genus_term}_{date}"
max_len = 5000  # Excluding too long (genomic) Sequences

additional_terms = []
struct = []

"""
if date_start == -1:
    date_start = 

if date_end != None:
    date_end = args.end
"""

# Make temporary directory if not exists

try:
    os.mkdir("./tmp")
except:
    pass


# Download all seqs from NCBI with genus Penicillium and save into json
if not (genus_term is None):
    Gen.NCBI_Download(email, genus_term, path_localgb)
elif not (accession_file is None):
    try:
        with open(accession_file, "r") as f:
            accession_list = [l.strip() for l in f.readlines()]
        Gen.NCBI_Downloadbyacclist(email, accession_list, path_localgb)
    except:
        Gen.Mes(f"Accession file {accession_file} is not valid!")
        raise Exception
else:
    Gen.Mes("Either genus name nor accession file is needed!")
    raise Exception


# Turn downloaded json to xlsx
Gen.jsontoxlsx(path_localgb, path_localgb_xlsx, max_len)

# Get seqs with term Korea
Gen.getseq(path_localgb, genus_term, path_out, additional_terms=additional_terms)
Gen.jsontransform(path_out + ".json", path_out + "_transformed.json")

Gen.jsontoxlsx(path_out + ".json", path_out + ".xlsx", max_len)
Gen.uni_jsontoxlsx(path_out + "_transformed.json", path_out + "_transformed.xlsx")

korean_acc = Gen.get_acc(path_out + "_transformed.json")

Gen.Mes(f"Total {len(korean_acc)} found")
classified_genes = Gen.classifier(path_out + "_transformed.json", path_out)
print(classified_genes)


# underlying module is for collecting nearby sequence to the genus, to prevent not collecing misidentified sequences
# however, blast is needed

"""
Gen.gene_seperator(classified_genes,path_out, email)

final_putative_list = [file for file in os.listdir(os.getcwd()) if file.startswith("Putative") and file.endswith(".json")]

Gen.json_merge(final_putative_list, f"Putative_{path_out}_All.json")
Gen.uni_jsontoxlsx(f"Putative_{path_out}_All.json", f"Putative_{path_out}_All.xlsx")
"""

# Build_DB(path_out+".fasta")
# BLASTn(path_out+".fasta",path_out+".fasta",10,"out.txt")
# BLAST_downloader(path_out+".fasta", path_out)


# type_file = SeqIO.parse(path_out, "fasta")
# for seq in type_file:
#   Struct_seq("TYPE_"+gene, seq, struct, "DB", "NaN")


# log_file.close()
