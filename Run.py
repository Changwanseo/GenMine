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


# For each genes, BLAST against NCBI DB, getoutput, reduction by local BLAST
for gene in classified_genes:
	if gene!="others" and gene!="genomic":
		
		print("Checkpoint 1")
		set_id = set()
		print(gene)
		list_foronlineblast = discrete_seqs(path_out+"_"+gene+".fasta")
		SeqIO.write(list_foronlineblast,f"onlineblastquery_{gene}.fasta","fasta")
		BLAST_downloader(f"onlineblastquery_{gene}.fasta", f"outonlineblast_{gene}.out")

		print("Checkpoint 2")
		for blastout in [file for file in os.listdir(os.getcwd()) if file.startswith(f"outonlineblast_{gene}.out")]:
			list_id = get_ids(blastout)
			set_id = set_id | set(list_id)

		print("Checkpoint 3")
		list_ID = list(set_id) 
		print(len(list_ID))

		print("Checkpoint 4")

		if len(list_ID) > 0:
			NCBI_Downloadbyacclist(email, list_ID, f"outonlineblast_{gene}.json")
			getseq_without(f"outonlineblast_{gene}.json",f"Byblast_{gene}", additional_terms=additional_terms, exceptional_terms=[genus_term])
			jsontransform(f"Byblast_{gene}.json",f"Byblast_transformed_{gene}.json")
			jsontoxlsx(f"Byblast_{gene}.json", f"Byblast_{gene}.xlsx", 3000)
			uni_jsontoxlsx(f"Byblast_transformed_{gene}.json", f"Byblast_transformed_{gene}.xlsx")

		else:
			Mes(f"No accessions available on {gene}")

for gene in classified_genes:
	non_detected_acc = []
	print(gene)
	if gene!="others" and gene!="genomic":
		acc_list = get_acc(f"Byblast_transformed_{gene}.json")
		for acc in acc_list:
			if acc in korean_acc:
				print(acc)
			else:
				non_detected_acc.append(acc)

		print(non_detected_acc)
		check = json_pick(f"Byblast_transformed_{gene}.json", non_detected_acc, f"Putative_{gene}.json")
		if check == True:
			uni_jsontoxlsx(f"Putative_{gene}.json", f"Putative_{gene}.xlsx")

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