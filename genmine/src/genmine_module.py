from Bio import Entrez
from Bio import SeqIO

# from Bio import pairwise2
# from Bio.pairwise2 import format_alignment

# from Bio.Blast.Applications import NcbiblastnCommandline
# from Bio.Blast import NCBIXML
from Bio.Seq import Seq
import time
from time import sleep
import pickle
import json
import re
from datetime import datetime
from itertools import repeat
import xmltodict
import os, sys, subprocess
import random
import re
import pickle
import pandas as pd
import shutil
import logging

# log_file = open("log.txt", "w")


# Logging Functions
def Time_now():
    now = datetime.now()
    Date_time = "%s-%s-%s %s:%s:%s " % (
        now.year,
        now.month,
        now.day,
        now.hour,
        now.minute,
        now.second,
    )
    return Date_time


# To both print to stdout and logfile
"""
def Mes(text):
    text = str(text)
    text = Time_now() + text
    log_file.write(text + "\n")
    print(text + "\n")
"""


def listtostring(list_in):
    if type(list_in) == type(["list"]):
        return ", ".join(list_in)
    else:
        return list_in


# Pretty print json
def print_json(x):
    print(json.dumps(x, indent=2))


def xml2dict(record):
    dict_type = xmltodict.parse(record)
    json_type = json.dumps(dict_type, indent=4)
    dict2_type = json.loads(json_type)
    return dict2_type


def delcomma(string):
    if string.endswith(","):
        return string[:-1]
    else:
        return string


# Main GenBank file downloadig function
def downloader(list_acc, path_tmp, out, cut=50):
    # Parse
    cnt = 0  # counter by 50
    cnt_all = len(list_acc)  # total number of records

    start_time = time.time()

    record_list = []

    if cnt_all == 0:
        logging.error("No accessions available in GenBank. Please check your input")
        return -1

    # Iterating with 50 items
    for i in range(int((len(list_acc) - 1) / cut) + 1):
        # Try to parse saved files
        if str(i) in [file for file in os.listdir(f"{path_tmp}")]:
            logging.info("Found saved")
            cnt += cut
            with open(f"{path_tmp}/{i}", "rb") as fp:
                data = pickle.load(fp)
                for record in data:
                    record_list.append(record)

        # for last last chunk
        elif i * cut + cut > len(list_acc):
            acc_string = ",".join(list_acc[i * cut :])
            sleep(0.3)
            cnt += len(list_acc[i * cut :])
            logging.info(
                f"{cnt}/{cnt_all} {round(100*cnt/cnt_all,2)}% {round(time.time()-start_time,2)}s"
            )

            # Parse xml
            # First trial
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=acc_string, rettype="gb", retmode="xml"
                )

            except:
                logging.warning("Server not responding, Requesting again...")

                # Second trial
                try:
                    sleep(10)
                    handle = Entrez.efetch(
                        db="nucleotide", id=acc_string, rettype="gb", retmode="xml"
                    )

                # Third trial
                except:
                    try:
                        sleep(100)
                        handle = Entrez.efetch(
                            db="nucleotide", id=acc_string, rettype="gb", retmode="xml"
                        )
                    # Final trial
                    except:
                        logging.warning(
                            "Last requesting due to connection error, will be take long (15 min)..."
                        )
                        sleep(900)
                        handle = Entrez.efetch(
                            db="nucleotide", id=acc_string, rettype="gb", retmode="xml"
                        )

            # Read parsed data
            pre_record = handle.read()
            json_record = xml2dict(pre_record)
            tmp_record_list = []

            # If only one result available, change record to list format
            if type(json_record["GBSet"]["GBSeq"]) is dict:
                json_record["GBSet"]["GBSeq"] = [json_record["GBSet"]["GBSeq"]]

            for record in json_record["GBSet"]["GBSeq"]:
                record_list.append(record)
                tmp_record_list.append(record)

            # Intermediate save for current data
            try:
                with open(f"{path_tmp}/{i}", "wb") as f:
                    pickle.dump(tmp_record_list, f)
            except:
                logging.error("[DEVELOPMENTAL ERROR] Data saving failed")
                raise Exception

        # non last chunk
        else:
            acc_string = ",".join(list_acc[i * cut : i * cut + cut])
            sleep(0.3)
            cnt += cut
            logging.info(
                f"{cnt}/{cnt_all} {round(100*cnt/cnt_all,2)}% {round(time.time()-start_time,2)}s"
            )

            # First trial
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=acc_string, rettype="gb", retmode="xml"
                )
            except:
                # Second trial
                logging.warning("Server not responding, Requesting again...")
                try:
                    sleep(10)
                    handle = Entrez.efetch(
                        db="nucleotide", id=acc_string, rettype="gb", retmode="xml"
                    )
                except:
                    # Third trial
                    try:
                        sleep(100)
                        handle = Entrez.efetch(
                            db="nucleotide", id=acc_string, rettype="gb", retmode="xml"
                        )
                    # Last trial
                    except:
                        logging.warning(
                            "Last requesting due to connection error, will be take long (15 min)..."
                        )
                        sleep(900)
                        handle = Entrez.efetch(
                            db="nucleotide", id=acc_string, rettype="gb", retmode="xml"
                        )

            # Read parsed data
            pre_record = handle.read()
            json_record = xml2dict(pre_record)
            tmp_record_list = []

            # If only one result available
            if type(json_record["GBSet"]["GBSeq"]) is dict:
                json_record["GBSet"]["GBSeq"] = [json_record["GBSet"]["GBSeq"]]

            for record in json_record["GBSet"]["GBSeq"]:
                record_list.append(record)
                tmp_record_list.append(record)

            try:
                with open(f"{path_tmp}/{i}", "wb") as f:
                    pickle.dump(tmp_record_list, f)
            except:
                logging.error("[DEVELOPMENTAL ERROR] Data saving failed")
                raise Exception

    with open(out, "w") as fp:
        json_term = json.dump(record_list, fp, indent=4)

    # If success, remove tmp files
    tmp_file_list = [file for file in os.listdir(f"{path_tmp}")]
    for file in tmp_file_list:
        os.remove(f"{path_tmp}/{file}")

    return 1


# Get accession list from GenBank
def ncbi_getacc(email, term, out):
    Entrez.email = email

    # Get all ID
    handle = Entrez.esearch(db="Nucleotide", term=term)
    record = Entrez.read(handle)

    sleep(5)
    handle = Entrez.esearch(
        db="Nucleotide", term=term, retmax=record["Count"], idtype="acc"
    )
    record = Entrez.read(handle)

    list_acc = record["IdList"]

    logging.info(f"Number of IDs: {len(list_acc)}")

    return list_acc


# To manage shotgun contig accessions
def entrez_query_generator(acc_list):
    # Shuffle acc_list to prevent errors in rerunning acc_list
    random.shuffle(acc_list)

    acc_tmp_list = []
    for acc in acc_list:
        # for normal genbank accession
        if re.fullmatch(
            r"(([A-Z]{1}[0-9]{5})(\.[0-9]{1}){0,1})|(([A-Z]{2}[\_]{0,1}[0-9]{6}){1}([\.][0-9]){0,1})",
            acc.strip(),
        ):
            acc_tmp_list.append(f"{acc.strip().split('.')[0]}[Accession]")
        # For shotgun sequences
        elif re.fullmatch(
            r"(([A-Z]{4}[0-9]{8})(\.[0-9]{1}){0,1})|(([A-Z]{6}[0-9]{9,})(\.[0-9]{1}){0,1})",
            acc.strip(),
        ):
            acc_tmp_list.append(acc.strip().split(".")[0] + "[Accession]")
        else:
            logging.error(
                f"[DEVELOPMENTAL ERROR] Bad accession {acc} entered entrez_query_generator. Please report your failed input to wan101010@snu.ac.kr"
            )
            raise Exception

    acc_string = " OR ".join(acc_tmp_list)

    return acc_string


def validate_accessions_binary_search(accession_list):
    """
    Validate GenBank accessions using binary search approach to efficiently
    identify invalid accessions.

    Args:
        accession_list: List of accession numbers to validate

    Returns:
        list: Valid accession numbers
    """

    # Binary search approach
    logging.info(f"Validating {len(accession_list)} accessions")

    try:
        sleep(0.3)
        # print(f"First of accession_list {accession_list[0]}")

        term = entrez_query_generator(accession_list)
        handle = Entrez.esearch(
            db="Nucleotide", term=term, retmax=2 * len(accession_list), idtype="acc"
        )

        http_response = str(handle.read())

        # Remove invalid ids
        if "ErrorList" in http_response:
            err_term = http_response.split("<ErrorList>")[1].split("</ErrorList>")[0]
            if "PhraseNotFound" in err_term:
                failed_acc = []
                tmp_terms = err_term.split("</PhraseNotFound>")

                for term in tmp_terms:
                    if "[Accession]" in term:
                        failed_acc.append(
                            term.split("<PhraseNotFound>")[1].split("[Accession]")[0]
                        )

            for acc in failed_acc:
                # logging.warning(f"GenBank accession {acc} is invalid")
                accession_list.remove(acc)

        return accession_list
    except:
        # print("Fail")
        if len(accession_list) == 1:
            logging.warning(f"GenBank accession {accession_list[0]} is invalid")
            return []
        else:
            # Split the list in half
            mid = len(accession_list) // 2
            left_half = accession_list[:mid]
            right_half = accession_list[mid:]

            # Recursively validate each half
            valid_left = validate_accessions_binary_search(left_half)
            valid_right = validate_accessions_binary_search(right_half)

            # Combine valid accessions from both halves
            return valid_left + valid_right


# filter accessions by regex
def filter_acc(acc_list, email, batch_size=100) -> list:
    """
    Filter and validate accessions with batching to prevent HTTP 500 errors.
    
    Args:
        acc_list: List of accession numbers
        email: Email for Entrez
        batch_size: Number of accessions to process per batch (default: 100)
    """
    # We can do this with pythonic expressions, but using iterative way for reporting
    return_acc_list = []
    return_shotgun_list = []
    
    # First, categorize accessions
    for acc in acc_list:
        if acc.strip() == "":
            pass
        # for normal genbank accession
        elif re.fullmatch(
            r"(([A-Z]{1}[0-9]{5})(\.[0-9]{1}){0,1})|(([A-Z]{2}[\_]{0,1}[0-9]{6}){1}([\.][0-9]){0,1})",
            acc.strip(),
        ):
            return_acc_list.append(acc.strip().split(".")[0])
        # For shotgun sequences
        elif re.fullmatch(
            r"(([A-Z]{4}[0-9]{8})(\.[0-9]{1}){0,1})|(([A-Z]{6}[0-9]{9,})(\.[0-9]{1}){0,1})",
            acc.strip(),
        ):
            return_shotgun_list.append(acc.strip().split(".")[0])
        else:
            logging.warning(
                f"Accession {acc} does not seems to be valid accession. Passing"
            )

    Entrez.email = email
    valid_acc_list = []

    if len(return_acc_list) > 0:
        logging.info(f"Processing {len(return_acc_list)} normal accessions in batches of {batch_size}...")
        
        # Process normal accessions in batches
        for i in range(0, len(return_acc_list), batch_size):
            batch = return_acc_list[i:i + batch_size]
            batch_num = (i // batch_size) + 1
            total_batches = (len(return_acc_list) + batch_size - 1) // batch_size
            
            logging.info(f"Processing batch {batch_num}/{total_batches} ({len(batch)} accessions)")
            
            try:
                query_term = entrez_query_generator(batch)
                logging.info(f"Query length: {len(query_term)} characters")
                
                handle = Entrez.esearch(
                    db="Nucleotide",
                    term=query_term,
                    retmax=2 * len(batch),
                    idtype="acc",
                )
                
                record = Entrez.read(handle)
                batch_results = [acc.split(".")[0] for acc in record["IdList"]]
                valid_acc_list.extend(batch_results)
                
                logging.info(f"Batch {batch_num} completed: found {len(batch_results)} valid accessions")
                
                # Check for missing accessions in this batch
                for acc in batch:
                    if acc not in batch_results:
                        logging.warning(f"GenBank accession {acc} cannot be found in GenBank. Please check misspelling or if the accessions are not opened yet.")
                
                # Be respectful to NCBI servers - add delay between batches
                sleep(0.5)
                
            except Exception as e:
                logging.warning(f"Error processing batch {batch_num}: {type(e).__name__}: {e}")
                # Continue with next batch instead of failing completely
                continue

    # Process shotgun accessions individually (as in original)
    if len(return_shotgun_list) > 0:
        logging.info(f"Processing {len(return_shotgun_list)} shotgun accessions individually...")
        for i, acc in enumerate(return_shotgun_list):
            try:
                handle = Entrez.esearch(
                    db="Nucleotide",
                    term=acc,
                    retmax=2,
                    idtype="acc",
                )
                record = Entrez.read(handle)
                shotgun_results = [_acc.split(".")[0] for _acc in record["IdList"]]
                valid_acc_list.extend(shotgun_results)
                
                if not acc in shotgun_results:
                    logging.warning(f"GenBank shotgun accession {acc} cannot be found in GenBank. Please check misspelling or if the accessions are not opened yet.")
                
                # Add small delay for shotgun processing too
                if i % 10 == 0:  # Every 10 accessions
                    sleep(0.1)
                    
            except Exception as e:
                logging.warning(f"Error processing shotgun accession {acc}: {e}")
                continue

    logging.info(f"Total valid accessions found: {len(valid_acc_list)}")
    return valid_acc_list


# GenBank record parser
# get wanted object from iterating dictionary
# if wanted object is in value of label
def retrieve(input_dict, obj, filter_list=[], add=" = ", default=""):
    obj_list = []
    for key in input_dict.keys():
        if key == obj:
            if not (type(input_dict[key]) is dict):
                return_list = []
                if type(input_dict[key]) in (str, int, float):
                    return_list.append(str(input_dict[key]))
                elif type(input_dict[key]) is list:
                    return_list += [str(i) for i in input_dict[key]]
                return return_list

        elif type(input_dict[key]) is dict:
            obj_list += retrieve(input_dict[key])

    out_list = []

    for o in obj_list:
        if not (o in filter_list):
            out_list.append(o)

    if len(out_list) == 0:
        return default
    else:
        return add.join(" = ")


# GenBank record parser
# get wanted object from iterating dictionary
# if wanted object and label is in parellel location
def retrieve_parallel(input_dict, label, label_value, obj):
    obj_list = []

    for key in input_dict.keys():
        if key == label and input_dict[label] == label_value:
            # print(input_dict)
            # find parallel keys
            for k in input_dict.keys():
                if k == obj:
                    if type(input_dict[obj]) in (str, int, float):
                        obj_list.append(str(input_dict[obj]))
                    elif type(input_dict[obj]) in (list, set):
                        obj_list += [str(i) for i in input_dict[obj]]

        elif type(input_dict[key]) is list:
            for o in input_dict[key]:
                if type(o) is dict:
                    obj_list += retrieve_parallel(o, label, label_value, obj)

        elif type(input_dict[key]) is dict:
            obj_list += retrieve_parallel(input_dict[key], label, label_value, obj)

    return obj_list


# change list output as string form
def format_list(input_list, filter_list=[], add=", ", default=""):
    out_list = []
    for o in input_list:
        if not (o in filter_list):
            out_list.append(o)

    if len(out_list) == 0:
        return default
    else:
        return add.join(out_list)


# json to xlsx for GenBank json file
def jsontoxlsx(json_in, xlsx, max_len=0):
    with open(json_in) as json_file:
        json_data = json.load(json_file)

    dict_all = {}
    for record in json_data:
        try:
            if (
                int(record["GBSeq_length"]) < max_len
            ):  # in order to get rid of genome data
                for key in record.keys():
                    if not (key in dict_all):
                        dict_all[key] = []
        except:
            # print_json(json_data)
            raise Exception

    for record in json_data:
        if int(record["GBSeq_length"]) < max_len:  # in order to get rid of genome data
            for key in dict_all.keys():
                try:
                    dict_all[key].append(str(record[key]))

                except:
                    dict_all[key].append("NaN")

    df = pd.DataFrame(dict_all)

    with pd.ExcelWriter(xlsx) as writer:
        df.to_excel(writer, index=False, sheet_name="Sheet 1")


# json to xlsx for universal json files
def uni_jsontoxlsx(json_in, xlsx):
    with open(json_in) as json_file:
        json_data = json.load(json_file)

    dict_all = {}
    for record in json_data:
        for key in record.keys():
            if not (key in dict_all):
                dict_all[key] = []

    for record in json_data:
        for key in dict_all.keys():
            try:
                dict_all[key].append(listtostring(record[key]))

            except:
                dict_all[key].append("NaN")

    df = pd.DataFrame(dict_all)

    with pd.ExcelWriter(xlsx) as writer:
        df.to_excel(writer, index=False, sheet_name="Sheet 1")


def jsontransform(json_in, out):  # transform to form easy to use
    # Get journal from given record json
    def Get_journal(record):
        state = 0
        journal = []
        department = []

        if "GBSeq_references" in record.keys():
            if "GBReference" in record["GBSeq_references"].keys():
                if type(record["GBSeq_references"]["GBReference"]) == dict:
                    if (
                        "submitted"
                        in record["GBSeq_references"]["GBReference"][
                            "GBReference_journal"
                        ].lower()
                    ):
                        department.append(
                            record["GBSeq_references"]["GBReference"][
                                "GBReference_journal"
                            ]
                        )
                    elif (
                        "unpublished"
                        in record["GBSeq_references"]["GBReference"][
                            "GBReference_journal"
                        ].lower()
                    ):
                        journal.append(
                            record["GBSeq_references"]["GBReference"][
                                "GBReference_journal"
                            ]
                        )
                    else:
                        journal.append(
                            record["GBSeq_references"]["GBReference"][
                                "GBReference_journal"
                            ]
                        )
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
                                    journal.remove(i)
                else:
                    raise Exception
            else:
                if "GBReferene_jounral" in record.keys():
                    if "submitted" in record["GBReference_jounal"].lower():
                        department.append(record["GBReference_jounal"])
                    else:
                        journal.append(record["GBReference_jounal"])
                else:
                    raise Exception

        else:
            return ["Unpublished"], []

        if journal == []:
            journal = ["Unpublished"]

        if department == []:
            pass

        return journal, department

    # Get paper title from given record json
    def Get_title(record):
        state = 0
        title = []
        if "GBSeq_references" in record:
            for reference in record["GBSeq_references"]["GBReference"]:
                if type(reference) == type({"dict": "dict"}):
                    try:
                        # print(len(reference))
                        if "GBReference_title" in reference:
                            if (
                                reference["GBReference_title"] == "Direct Submission"
                                and state <= 1
                            ):
                                state = 1  # only unpublished file exists
                                title.append("Direct Submission")

                            elif (
                                reference["GBReference_title"] == "Direct Submission"
                                and state >= 2
                            ):
                                pass

                            elif reference["GBReference_title"] != "Direct Submission":
                                state = 2
                                # print(title)
                                if "Direct Submission" in title:
                                    title.remove("Direct Submission")
                                title.append(reference["GBReference_title"])

                            else:  # unexpected input
                                logging.error(record["GBSeq_references"]["GBReference"])
                                logging.error(reference["GBReference_title"])
                                raise Exception
                        else:
                            if reference["GBReference_journal"] == "Unpublished":
                                state = 1  # only unpublished file exists
                                title.append("Direct Submission")
                            else:
                                # only journal record exists and no title exists
                                state = 2
                                title.append("Unknown")
                                # print(record)
                                # raise Exception

                    except:
                        logging.error(record)
                        raise Exception
        else:
            title = []

        return title

    # Get voucher from given record json
    def Get_voucher(record):
        obj_list = retrieve_parallel(
            input_dict=record,
            label="GBQualifier_name",
            label_value="specimen_voucher",
            obj="GBQualifier_value",
        )

        return format_list(input_list=obj_list, filter_list=[], add=" = ", default="")

    # Get type_material from given record json
    def Get_type_material(record):
        obj_list = retrieve_parallel(
            input_dict=record,
            label="GBQualifier_name",
            label_value="type_material",
            obj="GBQualifier_value",
        )

        return format_list(input_list=obj_list, filter_list=[], add=" = ", default="")

    # Get strain from given record json
    def Get_strain(record):
        obj_list = retrieve_parallel(
            input_dict=record,
            label="GBQualifier_name",
            label_value="strain",
            obj="GBQualifier_value",
        )

        return format_list(input_list=obj_list, filter_list=[], add=" = ", default="")

    # Get culture_collection from given record json
    def Get_culture_collection(record):
        obj_list = retrieve_parallel(
            input_dict=record,
            label="GBQualifier_name",
            label_value="culture_collection",
            obj="GBQualifier_value",
        )

        return format_list(input_list=obj_list, filter_list=[], add=" = ", default="")

    # Get note from given record json
    def Get_note(record):
        obj_list = retrieve_parallel(
            input_dict=record,
            label="GBQualifier_name",
            label_value="note",
            obj="GBQualifier_value",
        )

        return format_list(input_list=obj_list, filter_list=[], add=" = ", default="")

    # Get isolate from given record json
    def Get_isolate(record):
        obj_list = retrieve_parallel(
            input_dict=record,
            label="GBQualifier_name",
            label_value="isolate",
            obj="GBQualifier_value",
        )

        return format_list(input_list=obj_list, filter_list=[], add=" = ", default="")

    # Get clone from given record json
    def Get_clone(record):
        obj_list = retrieve_parallel(
            input_dict=record,
            label="GBQualifier_name",
            label_value="clone",
            obj="GBQualifier_value",
        )

        return format_list(input_list=obj_list, filter_list=[], add=" = ", default="")

    def Get_author(record):
        author_list = []

        if "GBSeq_references" in record.keys():
            if "GBReference" in record["GBSeq_references"].keys():
                if isinstance(record["GBSeq_references"]["GBReference"], dict):
                    if (
                        "GBReference_authors"
                        in record["GBSeq_references"]["GBReference"]
                    ):
                        if isinstance(
                            type(
                                record["GBSeq_references"]["GBReference"][
                                    "GBReference_authors"
                                ]["GBAuthor"]
                            ),
                            list,
                        ):
                            for author in record["GBSeq_references"]["GBReference"][
                                "GBReference_authors"
                            ]["GBAuthor"]:
                                if not (author in author_list):
                                    author_list.append(author)
                        elif (
                            type(
                                record["GBSeq_references"]["GBReference"][
                                    "GBReference_authors"
                                ]["GBAuthor"]
                            )
                            == str
                        ):
                            if not (
                                record["GBSeq_references"]["GBReference"][
                                    "GBReference_authors"
                                ]["GBAuthor"]
                                in author_list
                            ):
                                author_list.append(
                                    record["GBSeq_references"]["GBReference"][
                                        "GBReference_authors"
                                    ]["GBAuthor"]
                                )

                    elif (
                        "GBReference_consortium"
                        in record["GBSeq_references"]["GBReference"]
                    ):
                        if not (
                            record["GBSeq_references"]["GBReference"][
                                "GBReference_consortium"
                            ]
                            in author_list
                        ):
                            author_list.append(
                                record["GBSeq_references"]["GBReference"][
                                    "GBReference_consortium"
                                ]
                            )

                    else:
                        # print(record["GBSeq_references"]["GBReference"])
                        logging.warning("Failed to find authors")
                        return (
                            []
                        )  # in some of the records, no author exists. See NG_071242
                        # raise Exception

                elif isinstance(record["GBSeq_references"]["GBReference"], list):
                    # print("Multiple GBReference")
                    for reference in record["GBSeq_references"]["GBReference"]:
                        if "GBReference_authors" in reference:
                            if type(reference) == dict:
                                if "GBAuthor" in reference["GBReference_authors"]:
                                    # print(type(reference["GBReference_authors"]["GBAuthor"]))
                                    if (
                                        type(
                                            reference["GBReference_authors"]["GBAuthor"]
                                        )
                                        == list
                                    ):
                                        for author in reference["GBReference_authors"][
                                            "GBAuthor"
                                        ]:
                                            if not (author in author_list):
                                                author_list.append(author)
                                    elif (
                                        type(
                                            reference["GBReference_authors"]["GBAuthor"]
                                        )
                                        == str
                                    ):
                                        if (
                                            not (
                                                reference["GBReference_authors"][
                                                    "GBAuthor"
                                                ]
                                            )
                                            in author_list
                                        ):
                                            author_list.append(
                                                reference["GBReference_authors"][
                                                    "GBAuthor"
                                                ]
                                            )

                else:
                    logging.error(record)
                    raise Exception
            else:
                if "GBReferene_authors" in record.keys():
                    logging.error(record)
                    raise Exception
                else:
                    logging.error(record)
                    raise Exception

        else:
            logging.warning("No author information")

        author_list = list(set(author_list))

        return author_list

    with open(json_in) as json_file:
        json_data = json.load(json_file)

    json_temp = []

    # Transform to output format
    for record in json_data:
        dict_temp = {}

        if "GBSeq_sequence" in record or "GBSeq_feature-table" in record:
            # remove data without sequence
            dict_temp["acc"] = record["GBSeq_primary-accession"]
            dict_temp["length"] = record["GBSeq_length"]
            dict_temp["seqname"] = record["GBSeq_definition"]
            dict_temp["spname"] = record["GBSeq_organism"]
            dict_temp["uploader"] = Get_author(record)
            dict_temp["journal"], dict_temp["department"] = Get_journal(record)
            dict_temp["title"] = Get_title(record)
            dict_temp["upload_date"] = record["GBSeq_update-date"]  # GBSeq_create-date
            dict_temp["voucher"] = Get_voucher(record)
            dict_temp["type_material"] = Get_type_material(record)
            dict_temp["strain"] = Get_strain(record)
            dict_temp["culture_collection"] = Get_culture_collection(record)
            dict_temp["note"] = Get_note(record)
            dict_temp["isolate"] = Get_isolate(record)
            dict_temp["clone"] = Get_clone(record)

            dict_temp["primer"] = classification(record["GBSeq_definition"])
            json_temp.append(dict_temp)

            if "GBSeq_sequence" in record:  # normal sequence data
                dict_temp["seq"] = record["GBSeq_sequence"]
            elif "GBSeq_feature-table" in record:  # genomic data
                dict_temp["seq"] = "Genomic"

        else:
            logging.error(
                "[DEVELOPMETAL ERROR] Could not found GBSeq-sequence or GBSeq_feature-table from record. Please send your input to wan101010@snu.ac.kr"
            )
            print_json(record)
            raise Exception

    with open(out, "w") as fp:
        json_term = json.dump(json_temp, fp, indent=4)


def getseq(DB, out, additional_terms=[]):
    # terms are chosen in definition
    # additional terms are chosen in all parts

    temp_list = []
    with open(DB) as json_file:
        json_data = json.load(json_file)

    outfasta = open(out + ".fasta", "w")

    for record in json_data:
        if True:
            list_temp = []
            for key in record:
                if (
                    any(
                        term.lower() in str(record[key]).lower()
                        for term in additional_terms
                    )
                ) or len(additional_terms) == 0:
                    list_temp.append("1")
                else:
                    list_temp.append("0")
            if "1" in list_temp:
                if "GBSeq_definition" in record and "GBSeq_sequence" in record:
                    outfasta.write(">")
                    outfasta.write(record["GBSeq_definition"])
                    outfasta.write("\n")
                    outfasta.write(record["GBSeq_sequence"])
                    outfasta.write("\n")
                temp_list.append(record)

        else:
            pass
            # print(record)

    outjson = out + ".json"

    with open(outjson, "w") as fp:
        json_term = json.dump(temp_list, fp, indent=4)


def getseq_without(DB, out, additional_terms=[], exceptional_terms=[]):
    # terms are chosen in definition
    # additional terms are chosen in all parts

    temp_list = []
    with open(DB) as json_file:
        json_data = json.load(json_file)

    outfasta = open(out + ".fasta", "w")

    for record in json_data:
        list_temp = []
        for key in record:
            if type(record) != type("string"):
                if not (
                    any(
                        term.lower() in str(record[key]).lower()
                        for term in exceptional_terms
                    )
                ):
                    if any(
                        term.lower() in str(record[key]).lower()
                        for term in additional_terms
                    ):
                        # print(record[key])
                        list_temp.append("1")
                    else:
                        list_temp.append("0")
            else:
                pass
                # print(record)
        if "1" in list_temp:
            if "GBSeq_definition" in record and "GBSeq_sequence" in record:
                outfasta.write(">")
                outfasta.write(record["GBSeq_definition"])
                outfasta.write("\n")
                outfasta.write(record["GBSeq_sequence"])
                outfasta.write("\n")
            temp_list.append(record)

    outjson = out + ".json"

    with open(outjson, "w") as fp:
        json_term = json.dump(temp_list, fp, indent=4)


def seqrecordtodict(Seqrecord):
    dict_record = {
        "acc": "",
        "len": 0,
        "seqname": "",
        "spname": "",
        "uploader": [],
        "journal_things": [],
        "upload_date": "",
        "seq": "",
        "primer": "",
    }
    dict_record["acc"] = Seqrecord.id
    dict_record["len"] = len(Seqrecord)
    dict_record["seqname"] = Seqrecord.description

    dict_record["spname"] = Seqrecord.annotations["organism"]

    for Reference in Seqrecord.annotations["references"]:
        dict_record["journal_things"].append(Reference.title)
        dict_record["journal_things"].append(Reference.journal)

        for author in list(map(delcomma, Reference.authors.split(" "))):
            if not (author in dict_record["uploader"]):
                if author != "and":
                    dict_record["uploader"].append(author)

    dict_record["upload_date"] = Seqrecord.annotations["date"]
    dict_record["seq"] = str(Seqrecord.seq)
    dict_record["primer"] = classification(Seqrecord.description)
    return dict_record


def date_remover(string):
    try:
        return string.split(")")[1]
    except:
        # print(string)
        raise exception


def date_renamer(string):
    dict_mon = {
        "JAN": "01",
        "FEB": "02",
        "MAR": "03",
        "APR": "04",
        "MAY": "05",
        "JUN": "06",
        "JUL": "07",
        "AUG": "08",
        "SEP": "09",
        "OCT": "10",
        "NOV": "11",
        "DEC": "12",
    }
    spl = string.split("-")
    return f"{spl[2]}-{dict_mon[spl[1]]}-{spl[0]}"


def classification(description):
    if (
        "internal transcribed spacer 1" in description
        and "internal transcribed spacer 2" in description
    ):
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
    elif any(x in description.upper() for x in ["RPB2", "RNA POLYMERASE II", "RBP2"]):
        return "RPB2"
        print("RPB2")
    elif "actin" in description.lower() or "ACT" in description:
        return "ACT"
        print("ACT")
    elif "cox1" in description:
        return "cox1"
        print("cox1")
    elif "genomic sequence" in description:
        return "genomic"
        print("genomic")
    elif "COI" in description:
        return "COI(mt)"
        print("COI(mt)")
    elif (
        "28S ribosomal RNA" in description
        or "large subunit ribosomal RNA" in description
    ):
        return "LSU"
        print("LSU")
    elif (
        "small subunit ribosomal RNA" in description and "mitochondrial" in description
    ):
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
        # print(description)
        return "others"


def BLAST_downloader(fasta_in, blast_out):
    records = SeqIO.parse(fasta_in, "fasta")
    for i, record in enumerate(list(records)):
        # print(i)
        query = SeqIO.write(record, "temp.fasta", "fasta")
        os.system(
            f"blastn -out {blast_out}_{i}.xml -outfmt 5 -query temp.fasta -db /data/cwseo/BLAST_DB/public/nt -evalue 0.1 -num_threads 2"
        )
        # blastn_cline = ncbiblastnCommandline(query="temp.fasta",db="nr",evalue=0.1,outfmt=7,out=f"{blast_out}_{i}.xml")
        # print(blastn_cline)
        # blastn_cline()


def Build_DB(DB_fasta, out):
    cmd = "makeblastdb -in " + DB_fasta + " -out " + out + " -dbtype nucl"
    os.system(cmd)


def BLASTn(query, db, evalue, out):
    blastn_cline = ncbiblastnCommandline(
        query=query, db=db, evalue=evalue, outfmt=7, out=out
    )
    logging.info("BLAST: " + str(blastn_cline))
    blastn_cline()


def classifier(json_in, out):
    with open(json_in) as json_file:
        json_data = json.load(json_file)

    dict_fasta = {}

    for record in json_data:
        if not (record["primer"] in dict_fasta):
            dict_fasta[record["primer"]] = []
            dict_fasta[record["primer"]].append(f">{record['acc']}\n{record['seq']}\n")
        else:
            dict_fasta[record["primer"]].append(f">{record['acc']}\n{record['seq']}\n")

    new_keys = []

    for key in dict_fasta.keys():
        new_key = (
            key.replace(",", "_").replace(" ", "_").replace("(", "").replace(")", "")
        )
        new_keys.append(new_key)
        file = open(out + "_" + new_key + ".fasta", "w")
        for seq in dict_fasta[key]:
            file.write(seq)
        file.close()

    return new_keys


# get list of acc in transformed json
def get_acc(json_file):
    file = open(json_file, encoding="UTF-8")
    json_list = json.loads(file.read())

    list_acc = []

    for read in json_list:
        list_acc.append(read["acc"])

    return list_acc


# from json file, pick data with given acc and save as json
def json_pick(json_file, acc_list, out):
    file = open(json_file, encoding="UTF-8")
    json_list = json.loads(file.read())

    list_json = []

    for read in json_list:
        # print(read["acc"])
        if read["acc"] in acc_list:
            list_json.append(read)

    if len(list_json) == 0:
        logging.warning(
            f"All accessions in {json_file} has already found in genbank search"
        )
        return False

    else:
        with open(out, "w") as fp:
            json.dump(list_json, fp, indent=4)

        return True


def json_merge(json_list, out):
    all_records = []
    for json_file in json_list:
        file = open(json_file, encoding="UTF-8")
        tmp_list = json.loads(file.read())
        all_records += tmp_list

    with open(out, "w") as fp:
        json.dump(all_records, fp, indent=4)


# save download results to
def saver(path_work, name_out, max_len, additional_term):
    # Turn downloaded json to xlsx
    jsontoxlsx(f"{path_work}/{name_out}.json", f"{path_work}/{name_out}.xlsx", max_len)

    # Get sequence file with additional term
    getseq(
        DB=f"{path_work}/{name_out}.json",
        out=f"{path_work}/{name_out}",
        additional_terms=additional_term,
    )

    # Transform json into user friendly form
    jsontransform(
        f"{path_work}/{name_out}.json",
        f"{path_work}/{name_out}_transformed.json",
    )

    # Transform json to excel form
    jsontoxlsx(
        f"{path_work}/{name_out}.json",
        f"{path_work}/{name_out}.xlsx",
        max_len,
    )

    # Transform user friendly json to excel form
    uni_jsontoxlsx(
        f"{path_work}/{name_out}_transformed.json",
        f"{path_work}/{name_out}_transformed.xlsx",
    )

    # Parse finalized data
    term_acc = get_acc(f"{path_work}/{name_out}_transformed.json")

    logging.info(f"Total {len(term_acc)} found")
    classified_genes = classifier(
        f"{path_work}/{name_out}_transformed.json",
        f"{path_work}/{name_out}",
    )


# GenBank, download all by given term
# genus term is space seperated string of genus (eg. Penicillium Apiospora Alternaria )
def ncbi_download(
    email, genus_term, additional_term, name_out, path_work, path_tmp, max_len
):
    # select outgroup location
    path_localgb = f"{path_work}/{name_out}.json"
    path_localgb_xlsx = f"{path_work}/{name_out}.xlsx"

    # Setup email
    Entrez.email = email

    # Get all ID for given term
    handle = Entrez.esearch(db="Nucleotide", term=genus_term)
    record = Entrez.read(handle)

    # Get all GenBank records
    handle = Entrez.esearch(db="Nucleotide", term=genus_term, retmax=record["Count"])
    record = Entrez.read(handle)

    # Parse number of records
    list_acc = record["IdList"]

    logging.info(f"Number of IDs: {len(list_acc)}")

    # Download GenBank records
    status = downloader(
        list_acc=list_acc, path_tmp=path_tmp, out=f"{path_work}/{name_out}.json"
    )

    if not (status == -1):
        # save files
        saver(
            path_work=path_work,
            name_out=name_out,
            max_len=max_len,
            additional_term=additional_term,
        )


# Download by given list of accession
def ncbi_downloadbyacclist(email, list_acc, name_out, path_work, path_tmp, max_len):
    Entrez.email = email

    logging.info(f"Number of IDs: {len(list_acc)}")
    # Download GenBank records
    status = downloader(
        list_acc=list_acc, path_tmp=path_tmp, out=f"{path_work}/{name_out}.json"
    )

    if not (status == -1):
        # save files
        saver(
            path_work=path_work,
            name_out=name_out,
            max_len=max_len,
            additional_term=[],
        )
