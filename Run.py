from lib import GenMiner as Gen

import os

email = "wan101010@snu.ac.kr"
genus_term = "Umbelopsis" #Genus
date = "2020-04-22"
path_localgb = f"{genus_term}_{date}.json"
path_localgb_xlsx = f"{genus_term}_{date}.xlsx"
max_len = 5000 # Excluding too long (genomic) Sequences
additional_terms = []
struct = []
path_out = f"{genus_term}_{date}"



try:
	os.mkdir("./tmp")
except:
	pass
#Download all seqs from NCBI with genus Penicillium and save into json
Gen.NCBI_Download(email, genus_term, path_localgb)

# Turn downloaded json to xlsx
Gen.jsontoxlsx(path_localgb, path_localgb_xlsx, max_len)

# Get seqs with term Korea
Gen.getseq(path_localgb, genus_term, path_out, additional_terms=additional_terms)
Gen.jsontransform(path_out+".json",path_out+"_transformed.json")

Gen.jsontoxlsx(path_out+".json", path_out+".xlsx", 3000)
Gen.uni_jsontoxlsx(path_out+"_transformed.json",path_out+"_transformed.xlsx")

korean_acc = Gen.get_acc(path_out+"_transformed.json")

Gen.Mes(f"Total {len(korean_acc)} found")
classified_genes = Gen.classifier(path_out+"_transformed.json",path_out)
print(classified_genes)


Gen.gene_seperator(classified_genes,path_out, email)

final_putative_list = [file for file in os.listdir(os.getcwd()) if file.startswith("Putative") and file.endswith(".json")]

Gen.json_merge(final_putative_list, f"Putative_{genus_term}_All.json")
Gen.uni_jsontoxlsx(f"Putative_{genus_term}_All.json", f"Putative_{genus_term}_All.xlsx")

#Build_DB(path_out+".fasta")
#BLASTn(path_out+".fasta",path_out+".fasta",10,"out.txt")
#BLAST_downloader(path_out+".fasta", path_out)


#type_file = SeqIO.parse(path_out, "fasta")
#for seq in type_file:
#   Struct_seq("TYPE_"+gene, seq, struct, "DB", "NaN")


log_file.close()