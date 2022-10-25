def main():

    from . import GenMiner as Gen
    from lib.command import CommandParser
    from datetime import datetime, date

    import os

    args = CommandParser().get_args()

    date_start = -1
    date_end = -1

    if args.email != None:
        email = args.email
    else:
        print("No emails accepted, this might be bad for NCBI connection")

    genus_term = None
    accession_file = None

    if args.accession != None:
        accession_file = args.accession
    elif args.genus != None:
        genus_term = " ".join(args.genus)
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

    date = date.today().strftime("%Y-%m-%d")
    path_out = f"{genus_term}_{date}"
    path_localgb = f"{path_out}.json"
    path_localgb_xlsx = f"{path_out}.xlsx"

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
