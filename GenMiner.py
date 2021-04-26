from Bio import Entrez
from Bio import SeqIO
from Bio import pairwise2
from Bio.Blast import NCBIWWW
from Bio.pairwise2 import format_alignment
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
import time
from time import sleep
import pickle
from ete3 import Tree
import json
import re
from datetime import datetime
from itertools import repeat
import xmltodict
import os
import subprocess
import sys
import copy
import random
import pickle
import pandas as pd
import shutil
import multiprocessing as mp


log_file = open("log.txt","w")

#Logging Functions
def Time_now():
	now = datetime.now()
	Date_time = ("%s-%s-%s %s:%s:%s " % (now.year,now.month,now.day,now.hour,now.minute,now.second))
	return Date_time

def Mes(text):
	text = str(text)
	text = Time_now()+text
	log_file.write(text+"\n")
	print(text+"\n")

def listtostring(list_in):
	if type(list_in) == type(["list"]):
		return ", ".join(list_in)
	else:
		return list_in

def NCBI_Download(Email, term, out):

	def xml2dict(record):
		dict_type = xmltodict.parse(record)
		json_type = json.dumps(dict_type, indent=4)
		dict2_type= json.loads(json_type)
		return dict2_type

	Entrez.email = Email

	#Get all ID
	handle = Entrez.esearch(db="Nucleotide", term=term)
	record = Entrez.read(handle)
	#Mes(str(record))

	handle =  Entrez.esearch(db="Nucleotide", term=term, retmax=record['Count'])
	record = Entrez.read(handle)

	list_ID = record['IdList']

	print(len(list_ID))

	Mes(f"Number of IDs: {len(list_ID)}")

	#Parse
	cnt = 0
	cut = 50 #parse by 50 sequences
	cnt_all = len(list_ID)

	start_time = time.time()

	record_list = []

	#print(int((len(list_ID)-1)/cut)+1)

	for i in range(int((len(list_ID)-1)/cut)+1):
		if str(i) in [file for file in os.listdir("./tmp")]:
			Mes("Found saved")
			cnt+=cut
			with open(f"./tmp/{i}", "rb") as fp:
				data = pickle.load(fp)
				for record in data:
					record_list.append(record)

		elif i*cut+cut > len(list_ID):
			ID_string = ",".join(list_ID[i*cut:])
			sleep(0.3)
			cnt+=len(list_ID[i*cut:])
			Mes(f'{cnt}/{cnt_all} {100*cnt/cnt_all}% {time.time()-start_time}s')

			try: 
				handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
			except:
				Mes("Requesting again...")
				try:
					sleep(10)
					handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
				except:
					try:
						sleep(100)
						handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
					except:
						sleep(1000)
						handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")

			pre_record = handle.read()
			json_record = xml2dict(pre_record)
			tmp_record_list = []

			for record in json_record['GBSet']['GBSeq']:
				record_list.append(record)
				tmp_record_list.append(record)

			try:
				with open(f"./tmp/{i}", "wb") as f:
					pickle.dump(tmp_record_list,f)
			except:
				Mes("Saving Error")
				raise Exception



		else:
			Mes(i)
			ID_string = ",".join(list_ID[i*cut:i*cut+cut])
			sleep(0.3)
			cnt+=cut
			Mes(f'{cnt}/{cnt_all} {100*cnt/cnt_all}% {time.time()-start_time}s')

			try: 
				handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
			except:
				Mes("Requesting again...")
				try:
					sleep(10)
					handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
				except:
					try:
						sleep(100)
						handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
					except:
						sleep(1000)
						handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")

			pre_record = handle.read()
			json_record = xml2dict(pre_record)
			tmp_record_list = []

			for record in json_record['GBSet']['GBSeq']:
				record_list.append(record)
				tmp_record_list.append(record)

			try:
				with open(f"./tmp/{i}", "wb") as f:
					pickle.dump(tmp_record_list,f)
			except:
				Mes("Saving Error")
				raise Exception

	with open(out, "w") as fp:
		json_term = json.dump(record_list, fp, indent= 4)

def NCBI_Downloadbyacclist(Email, list_ID, out):

	#print(list_ID)

	def xml2dict(record):
		dict_type = xmltodict.parse(record)
		json_type = json.dumps(dict_type, indent=4)
		dict2_type= json.loads(json_type)
		return dict2_type

	Entrez.email = Email

	Mes(f"Number of IDs: {len(list_ID)}")

	#Parse
	cnt = 0
	cut = 50 #parse by 50 sequences
	cnt_all = len(list_ID)

	start_time = time.time()

	record_list = []

	for i in range(int((len(list_ID)-1)/cut)+1):
		if i*cut+cut > len(list_ID):
			Mes(i)
			ID_string = ",".join(list_ID[i*cut:])
			sleep(0.3)
			cnt+=cut
			Mes(f'{cnt}/{cnt_all} {100*cnt/cnt_all}% {time.time()-start_time}s')

			try: 
				handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
			except: #Retry
				Mes("Requesting again...")
				try:
					sleep(10)
					handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
				except:
					sleep(100)
					handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")

			pre_record = handle.read()
			json_record = xml2dict(pre_record)
			for record in json_record['GBSet']['GBSeq']:
				record_list.append(record)

			try:
				with open(f"./tmp/{i}", "wb") as f:
					pickle.dump(record,f)
			except:
				Mes("Saving Error")
				raise Exception

		else: 
			Mes(i)
			ID_string = ",".join(list_ID[i*cut:i*cut+cut])
			sleep(0.3)
			cnt+=cut
			Mes(f'{cnt}/{cnt_all} {100*cnt/cnt_all}% {time.time()-start_time}s')

			try: 
				handle = Entrez.efetch(db="nucleotide",id=ID_string,rettype="gb",retmode="xml")
			except: #Retry
				Mes("Requesting again...")
				try:
					sleep(10)
					handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")
				except:
					sleep(100)
					handle = Entrez.efetch(db="nucleotide",id=ID_string, rettype="gb",retmode="xml")

			pre_record = handle.read()
			json_record = xml2dict(pre_record)
			
			for record in json_record['GBSet']['GBSeq']:
				record_list.append(record)
				print(len(record_list))
			
			try:
				with open(f"./tmp/{i}", "wb") as f:
					pickle.dump(record,f)
			except:
				Mes("Saving Error")
				raise Exception

	with open(out, "w") as fp:
		json_term = json.dump(record_list, fp, indent= 4)

def jsontoxlsx(json_in, xlsx, max_len=0):
	with open(json_in) as json_file:
		json_data = json.load(json_file)

	dict_all = {}
	for record in json_data:
		if int(record["GBSeq_length"]) < max_len: #in order to get rid of genome data
			for key in record.keys():
				if not(key in dict_all):
					dict_all[key] = []

	for record in json_data:
		if int(record["GBSeq_length"]) < max_len: #in order to get rid of genome data
			for key in dict_all.keys():
				try:
					dict_all[key].append(str(record[key]))

				except:
					dict_all[key].append("NaN")

	df = pd.DataFrame(dict_all)

	with pd.ExcelWriter(xlsx) as writer:
		df.to_excel(writer, sheet_name = "Sheet 1")

def uni_jsontoxlsx(json_in, xlsx):
	with open(json_in) as json_file:
		json_data = json.load(json_file)

	dict_all = {}
	for record in json_data:
		for key in record.keys():
			if not(key in dict_all):
				dict_all[key] = []

	for record in json_data:
		for key in dict_all.keys():
			try:
				dict_all[key].append(listtostring(record[key]))

			except:
				dict_all[key].append("NaN")

	df = pd.DataFrame(dict_all)

	with pd.ExcelWriter(xlsx) as writer:
		df.to_excel(writer, sheet_name = "Sheet 1")

def jsontransform(json_in, out):

	def Get_journal(record):

		print("----------------------------")

		state = 0 
		journal = []
		department = []

		if "GBSeq_references" in record.keys():
			if "GBReference" in record["GBSeq_references"].keys():
				if type(record["GBSeq_references"]["GBReference"]) == dict:
					if "submitted" in record["GBSeq_references"]["GBReference"]["GBReference_journal"].lower():
						department.append(record["GBSeq_references"]["GBReference"]["GBReference_journal"])
					elif "unpublished" in record["GBSeq_references"]["GBReference"]["GBReference_journal"].lower():
						journal.append(record["GBSeq_references"]["GBReference"]["GBReference_journal"])
					else:
						journal.append(record["GBSeq_references"]["GBReference"]["GBReference_journal"])
						for i in journal:
							if "unpublished" in i.lower():
								jounral.remove(i)

				elif type(record["GBSeq_references"]["GBReference"]) == list:
					for reference in record["GBSeq_references"]["GBReference"]:
						if "submitted" in reference["GBReference_journal"].lower():
							department.append(reference["GBReference_journal"])
						elif "unpublished" in reference["GBReference_journal"].lower():
							journal.append(reference["GBReference_journal"])
						else:
							journal.append(reference["GBReference_journal"])
							for i in journal:
								if "unpublished" in i.lower():
									jounral.remove(i)
					
				else:
					print(record)
					raise Exception
				#for reference in record["GBSeq_references"]["GBReference"]:
			else:
				if "GBReferene_jounral" in record.keys():
					if "submitted" in record["GBReference_jounal"].lower():
						department.append(record["GBReference_jounal"])
					else:
						journal.append(record["GBReference_jounal"])
				else:
					print(record)
					raise Exception

		else:
			print(record)
			raise Exception

		

		if journal == []:
			journal = ["Unpublished"]

		if department == []:
			pass
			#print(json.dumps(record["GBSeq_references"], indent=2))
			#raise Exception


		print(f"journal: {journal}")
		print(f"department: {department}")


		return journal, department

	def Get_title(record):

		state = 0 
		title = []

		for reference in record["GBSeq_references"]["GBReference"]:
			if type(reference) == type({"dict":"dict"}):
				if reference["GBReference_title"] == "Direct Submission" and state <= 1:
					state = 1
					title.append("Direct Submission")

				elif reference["GBReference_title"] == "Direct Submission" and state >= 2:
					pass

				elif reference["GBReference_title"] != "Direct Submission":
					state = 2
					#print(title)
					if "Direct Submission" in title:
						title.remove("Direct Submission")
					title.append(reference["GBReference_title"])

				else: #unexpected input
					print(record["GBSeq_references"]["GBReference"])
					print(reference["GBReference_title"])
					raise Exception

		return title

	def Get_author(record):

		author_list = []

		if "GBSeq_references" in record.keys():
			print("Found GBSeq_references")
			if "GBReference" in record["GBSeq_references"].keys():
				print("Found GBReference")
				if type(record["GBSeq_references"]["GBReference"]) == dict:
					print("Single GBReference")
					if "GBReference_authors" in record["GBSeq_references"]["GBReference"]:
						if type(record["GBSeq_references"]["GBReference"]["GBReference_authors"]["GBAuthor"]) == list:
							for author in record["GBSeq_references"]["GBReference"]["GBReference_authors"]["GBAuthor"]:
								if not(author in author_list):
									author_list.append(author)
						elif type(record["GBSeq_references"]["GBReference"]["GBReference_authors"]["GBAuthor"]) == str:
							if not(record["GBSeq_references"]["GBReference"]["GBReference_authors"]["GBAuthor"] in author_list):
								author_list.append(record["GBSeq_references"]["GBReference"]["GBReference_authors"]["GBAuthor"])

					elif "GBReference_consortium" in record["GBSeq_references"]["GBReference"]:
						if not(record["GBSeq_references"]["GBReference"]["GBReference_consortium"] in author_list):
							author_list.append(record["GBSeq_references"]["GBReference"]["GBReference_consortium"])
						
					else:
						print(record["GBSeq_references"]["GBReference"])
						print("Failed to find authors")
						raise Exception

				elif type(record["GBSeq_references"]["GBReference"]) == list:
					print("Multiple GBReference")
					for reference in record["GBSeq_references"]["GBReference"]:
						if "GBReference_authors" in reference:
							if type(reference) == dict:
								if "GBAuthor" in reference["GBReference_authors"]:
									#print(type(reference["GBReference_authors"]["GBAuthor"]))
									if type(reference["GBReference_authors"]["GBAuthor"]) == list:
										for author in reference["GBReference_authors"]["GBAuthor"]:
											if not(author in author_list):
												author_list.append(author)
									elif type(reference["GBReference_authors"]["GBAuthor"]) == str:
										if not(reference["GBReference_authors"]["GBAuthor"]) in author_list:
											author_list.append(reference["GBReference_authors"]["GBAuthor"])	
					
				else:
					print(record)
					raise Exception
			else:
				if "GBReferene_authors" in record.keys():
					print(record)
					raise Exception
				else:
					print(record)
					raise Exception

		else:
			print(record)
			raise Exception


		for reference in record["GBSeq_references"]["GBReference"]:
			#print(json.dumps(reference, indent=2))
			if "GBReference_authors" in reference:
				if type(reference) == type({"dict":"dict"}):
					if "GBAuthor" in reference["GBReference_authors"]:
						#print(type(reference["GBReference_authors"]["GBAuthor"]))
						if type(reference["GBReference_authors"]["GBAuthor"]) == type(["l","i","s","t"]):
							for author in reference["GBReference_authors"]["GBAuthor"]:
								if not(author in author_list):
									author_list.append(author)
						elif type(reference["GBReference_authors"]["GBAuthor"]) == type("string"):
							author_list.append(reference["GBReference_authors"]["GBAuthor"])

		author_list = list(set(author_list))
		print(f"author: {author_list}")

		if author_list == []:
			print(json.dumps(record, indent = 2))
			raise Exception

		return author_list

	with open(json_in) as json_file:
		json_data = json.load(json_file)

	json_temp = []

	for record in json_data:
		if "GBSeq_sequence" in record: #remove data without sequence
			#print(json.dumps(record, indent=2))
			#print("---------------------------------")
			dict_temp = {"acc":"","length":"", "seqname":"", "spname":"", "uploader":[],"journal":[], "department":[],"title":"","upload_date":"", "seq":""}
			dict_temp["acc"] = record["GBSeq_locus"]
			dict_temp["length"] = record["GBSeq_length"]
			dict_temp["seqname"] = record["GBSeq_definition"]
			dict_temp["spname"] = record["GBSeq_organism"]
			dict_temp["uploader"] = Get_author(record)
			dict_temp["journal"], dict_temp["department"] = Get_journal(record)
			dict_temp["title"] = Get_title(record)
			dict_temp["upload_date"] =  record["GBSeq_update-date"] # GBSeq_create-date
			dict_temp["seq"] = record["GBSeq_sequence"]
			dict_temp["primer"] = classification(record["GBSeq_definition"])
			json_temp.append(dict_temp)

	with open(out, "w") as fp:
		json_term = json.dump(json_temp, fp, indent= 4)

def getseq(DB, gene, out, additional_terms = []):

	#terms are chosen in definition
	#additional terms are chosen in all parts

	temp_list = [] 
	with open(DB) as json_file:
		json_data = json.load(json_file)

	outfasta = open(out+".fasta","w")

	for record in json_data:
		if any(gene.lower() in str(record[term]).lower() for term in ["GBSeq_definition", "GBSeq_source", "GBSeq_organism", "GBSeq_taxonomy"]):
			list_temp = []
			for key in record:
				#print(str(record[key]).lower())
				if (any(term.lower() in str(record[key]).lower() for term in additional_terms)):
					#print(record[key])
					list_temp.append("1")
				else:
					list_temp.append("0")
			if ("1" in list_temp):
				if "GBSeq_definition" in record and "GBSeq_sequence" in record:
					outfasta.write(">")
					outfasta.write(record["GBSeq_definition"])
					outfasta.write("\n")
					outfasta.write(record["GBSeq_sequence"])
					outfasta.write("\n")
				temp_list.append(record)

	outjson = out+".json"

	with open(outjson, "w") as fp:
		json_term = json.dump(temp_list, fp, indent= 4)

def getseq_without(DB, out, additional_terms = [], exceptional_terms = []):

	#terms are chosen in definition
	#additional terms are chosen in all parts

	temp_list = [] 
	with open(DB) as json_file:
		json_data = json.load(json_file)

	outfasta = open(out+".fasta","w")

	for record in json_data:
		list_temp = []
		for key in record:
			if type(record) != type("string"):
				if not(any(term.lower() in str(record[key]).lower() for term in exceptional_terms)):
					if (any(term.lower() in str(record[key]).lower() for term in additional_terms)):
						#print(record[key])
						list_temp.append("1")
					else:
						list_temp.append("0")
			else:
				print(record)
		if ("1" in list_temp):
			if "GBSeq_definition" in record and "GBSeq_sequence" in record:
				outfasta.write(">")
				outfasta.write(record["GBSeq_definition"])
				outfasta.write("\n")
				outfasta.write(record["GBSeq_sequence"])
				outfasta.write("\n")
			temp_list.append(record)

	outjson = out+".json"

	with open(outjson, "w") as fp:
		json_term = json.dump(temp_list, fp, indent= 4)

def delcomma(string):
	if string.endswith(","):
		return string[:-1]
	else:
		return string

def seqrecordtodict(Seqrecord):

	dict_record = {"acc":"","len":0, "seqname":"", "spname":"", "uploader":[],"journal_things":[], "upload_date":"", "seq":"", "primer":""}
	dict_record["acc"] = Seqrecord.id
	dict_record["len"] = len(Seqrecord)
	dict_record["seqname"] = Seqrecord.description

	dict_record["spname"] = Seqrecord.annotations["organism"]

	for Reference in Seqrecord.annotations["references"]:

		dict_record["journal_things"].append(Reference.title)
		dict_record["journal_things"].append(Reference.journal)

		for author in list(map(delcomma,Reference.authors.split(" "))):
			if not(author in dict_record["uploader"]):
				if author!="and":
					dict_record["uploader"].append(author)

	dict_record["upload_date"] = Seqrecord.annotations["date"]
	dict_record["seq"] = str(Seqrecord.seq)
	dict_record["primer"] = classification(Seqrecord.description)
	return dict_record

def date_remover(string):
	try:
		return string.split(")")[1]
	except:
		print(string)
		raise exception
	
def date_renamer(string):
	dict_mon = {"JAN":"01","FEB":"02","MAR":"03","APR":"04","MAY":"05","JUN":"06","JUL":"07","AUG":"08","SEP":"09","OCT":"10","NOV":"11","DEC":"12"}
	spl = string.split("-")
	return f"{spl[2]}-{dict_mon[spl[1]]}-{spl[0]}"


def classification(description):

	if "internal transcribed spacer 1" in description and "internal transcribed spacer 2" in description:
		return "ITS1, ITS2"
		print("ITS1, ITS2")
	elif "internal transcribed spacer 1" in description:
		return "ITS1"
		print("ITS1")
	elif "internal transcribed spacer 2" in description:
		return "ITS2"
		print("ITS2") 
	elif "ITS" in description:
		return "ITS"
		print("ITS")
	elif any(x in description for x in ["beta-tubulin", "beta tubulin", "beta-tublin"]):
		return "TUB"
		print("TUB")
	elif "calmodulin" in description:
		return "CMD"
		print("CMD")
	elif "RPB2" in description.upper():
		return "RPB2"
		print("RPB2") 
	elif "cox1" in description:
		return "cox1"
		print("cox1") 
	elif "genomic sequence" in description:
		return "genomic"
		print("genomic") 
	elif "COI" in description:
		return "COI(mt)"
		print("COI(mt)") 
	elif "28S ribosomal RNA" in description or "large subunit ribosomal RNA" in description:
		return "LSU"
		print("LSU") 
	elif "small subunit ribosomal RNA" in description and "mitochondrial" in description:
		return "SSU(mt)"
		print("SSU(mt)")
	elif "18S ribosomal RNA" in description:
		return "SSU"
		print("SSU") 
	elif "elongation factor" in description:
		return "EF"
		print("EF") 
	elif "xynY" in description:
		return "xynY"
		print("xynY")
	elif "abfB" in description:
		return "abfB"
		print("abfB")
	elif "CHS3" in description:
		return "CHS3"
		print("CHS3")
	elif "CHS2" in description:
		return "CHS2"
		print("CHS2")
	elif "CHS1" in description:
		return "CHS1"
		print("CHS1")   
	elif "CHS4" in description:
		return "CHS4"
		print("CHS4")
	elif "PJ3" in description:
		return "PJ3"
		print("PJ3")
	elif "mitochond" in description:
		return "mt_others"
		print("mt_others")
	else:
		print(description)
		return "others"

def BLAST_downloader(fasta_in, blast_out):

	records = SeqIO.parse(fasta_in,"fasta")
	for i,record in enumerate(list(records)):
		print(i)
		#result_handle = NCBIWWW.qblast("blastn","nt",record.format("fasta"), hitlist_size=100000000,perc_ident=80)
		query = SeqIO.write(record,"temp.fasta","fasta")
		os.system(f"blastn -out {blast_out}_{i}.xml -outfmt 5 -query temp.fasta -db /data/cwseo/BLAST_DB/public/nt -evalue 0.1 -num_threads 2")
		#blastn_cline = NcbiblastnCommandline(query="temp.fasta",db="nr",evalue=0.1,outfmt=7,out=f"{blast_out}_{i}.xml")
		#print(blastn_cline)
		#blastn_cline()

def Build_DB(DB_fasta, out):

	cmd = 'makeblastdb -in '+ DB_fasta +' -out '+ out +" -dbtype nucl"
	os.system(cmd)

def BLASTn(query,db,evalue,out):
	blastn_cline = NcbiblastnCommandline(query=query,db=db,evalue=evalue,outfmt=7,out=out)
	Mes("BLAST: "+str(blastn_cline))
	blastn_cline()

def classifier(json_in,out):

	with open(json_in) as json_file:
		json_data = json.load(json_file)

	dict_fasta = {}

	for record in json_data:
		if not(record["primer"] in dict_fasta):
			dict_fasta[record["primer"]] = []
			dict_fasta[record["primer"]].append(f">{record['acc']}\n{record['seq']}\n")
		else:
			dict_fasta[record["primer"]].append(f">{record['acc']}\n{record['seq']}\n")

	new_keys = []

	for key in dict_fasta.keys():
		new_key = key.replace(",","_").replace(" ","_").replace("(","").replace(")","")
		new_keys.append(new_key)
		file = open(out+"_"+new_key+".fasta","w")
		for seq in dict_fasta[key]:
			file.write(seq)
		file.close()

	return new_keys

def discrete_seqs(fasta_file):

	Build_DB(fasta_file, "temp_db")

	fasta_list = list(SeqIO.parse(fasta_file, "fasta"))

	#print(fasta_list)

	SeqIO.write(fasta_list[0],"temp_query.fasta","fasta")
	BLASTn("temp_query.fasta", "temp_db", evalue=10, out="temp_blastout")
	
	#query
	list_id_out = [fasta_list[0]]

	#leftovers for next iteration
	list_id = get_ids2("temp_blastout")
	
	temp_fasta_list = []

	if len(list_id) > 2:

		#Save selected fasta to temp
		for fasta in fasta_list:
			for fasta_id in list_id:
				if fasta_id in fasta.description:
					temp_fasta_list.append(fasta)

		SeqIO.write(temp_fasta_list, "temp_db.fasta", "fasta")
		return list_id_out + discrete_seqs("temp_db.fasta")

	else:
		return list_id_out


def get_blast_result(blast_result):
	pass

#for local BLAST result
def get_ids2(blast_result):

	cutoff = 10

	list_id = []

	file = open(blast_result,"r")
	lines = file.readlines()

	for line in lines:
		if not(line.startswith("#")):
			if float(line.split("\t")[2]) < 100 - cutoff:
				#print(line)
				list_id.append(str(line.split("\t")[1]))

	return list_id

#for web BLAST result
def get_ids(file_xml):


	file = open(file_xml,"r")
	lines = file.readlines()

	list_id_num = []
	list_id = []

	for num, line in enumerate(lines):
		if "<Hit_id>" in line:
			list_id_num.append(num)

	for i in range(len(list_id_num)):
		try:
			flag = 0 
			for line in lines[list_id_num[i]:list_id_num[i+1]]:
				if "Penicillium" in line or "penicillium" in line:
					flag = 1

			if flag == 0:
				list_id.append(lines[list_id_num[i]].split("<Hit_id>")[1].split("</Hit_id>")[0].split("|")[1])

		except:
			flag = 0 
			for line in lines[list_id_num[i]:]:
				if "Penicillium" in line or "penicillium" in line:
					flag = 1

			if flag == 0:
				list_id.append(lines[list_id_num[i]].split("<Hit_id>")[1].split("</Hit_id>")[0].split("|")[1])

	for line in lines:
		if "<Hit_id>" in line:
			list_id.append(line.split("<Hit_id>")[1].split("</Hit_id>")[0].split("|")[1])

	file.close()
	print(len(list_id))
	return list_id

#get list of acc in transformed json
def get_acc(json_file):

	file = open(json_file,encoding="UTF-8")
	json_list = json.loads(file.read())

	list_acc = []

	for read in json_list:
		list_acc.append(read["acc"])

	return list_acc

#from json file, pick data with given acc and save as json
def json_pick(json_file, acc_list, out):

	file = open(json_file, encoding="UTF-8")
	json_list = json.loads(file.read())

	list_json = []

	for read in json_list:
		print(read["acc"])
		if read["acc"] in acc_list:
			list_json.append(read)

	if len(list_json) == 0:
		Mes(f"All accessions in {json_file} has already found in genbank search")
		return False

	else:
		with open(out, "w") as fp:
			json.dump(list_json, fp, indent= 4)

		return True

def json_merge(json_list, out):

	all_records = []
	for json_file in json_list:
		file = open(json_file, encoding="UTF-8")
		tmp_list = json.loads(file.read())
		all_records += tmp_list

	with open(out, "w") as fp:
		json.dump(all_records, fp, indent = 4)


email = "wan101010@snu.ac.kr"
genus_term = "Penicillium" #Genus`
path_localgb = "Penicillium_2020-04-22.json"
path_localgb_xlsx = "Penicillium_2020-04-22.xlsx"
max_len = 5000 # Excluding too long (genomic) Sequences
additional_terms = ["korea"]
struct = []
path_out = "Penicililum_Korea_2020-04-22"



try:
	os.mkdir("./tmp")
except:
	pass
#Download all seqs from NCBI with genus Penicillium and save into json
NCBI_Download(email, genus_term, path_localgb)

# Turn downloaded json to xlsx
jsontoxlsx(path_localgb, path_localgb_xlsx, max_len)

# Get seqs with term Korea
getseq(path_localgb, genus_term, path_out, additional_terms=additional_terms)
jsontransform(path_out+".json",path_out+"_transformed.json")

jsontoxlsx(path_out+".json", path_out+".xlsx", 3000)
uni_jsontoxlsx(path_out+"_transformed.json",path_out+"_transformed.xlsx")

korean_acc = get_acc(path_out+"_transformed.json")

Mes(f"Total {len(korean_acc)} found")
classified_genes = classifier(path_out+"_transformed.json",path_out)
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

json_merge(final_putative_list, f"Putative_{genus_term}_All.json")
uni_jsontoxlsx(f"Putative_{genus_term}_All.json", f"Putative_{genus_term}_All.xlsx")

#Build_DB(path_out+".fasta")
#BLASTn(path_out+".fasta",path_out+".fasta",10,"out.txt")
#BLAST_downloader(path_out+".fasta", path_out)


#type_file = SeqIO.parse(path_out, "fasta")
#for seq in type_file:
#   Struct_seq("TYPE_"+gene, seq, struct, "DB", "NaN")


log_file.close()