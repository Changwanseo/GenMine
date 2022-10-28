def main():

    from . import GenMine as Gen
    from lib.command import CommandParser
    from datetime import datetime, date

    import os

    args = CommandParser().get_args()

    # Managing arguments

    # default setup
    date_start = -1
    date_end = -1
    genus_term = None
    accession_file = None
    max_len = 5000  # Excluding too long (genomic) Sequences
    additional_term = []
    struct = []

    # Change setups by input
    if args.email != None:
        email = args.email
    else:
        email = None
        print("No emails accepted, this might be bad for NCBI connection")

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

    # make output directory name
    if not (genus_term is None):
        name_out = f"{genus_term}"
    elif not (accession_file is None):
        name_out = accession_file.replace("\\", "/").split("/")[-1].split(".")[0]
    else:
        Gen.Mes("Either genus name nor accession file is needed!")
        raise Exception

    # working directory
    path_work = f"{os.getcwd()}/{name_out}"

    # Check with previous run
    if os.path.isdir(path_work):
        print("Found previous run. Continue")
    else:
        try:
            os.mkdir(path_work)
        except:
            Gen.Mes(
                f"Could not create output location {path_work}. Check permissions or if previous run file exists"
            )
            raise Exception

    # Date options, should be edited
    """
    if date_start == -1:
        date_start = 

    if date_end != None:
        date_end = args.end
    """

    # Make temporary directory if not exists
    path_tmp = f"{path_work}/tmp"
    try:
        os.mkdir(path_tmp)
    except:
        pass

    # Download all seqs from NCBI with genus Penicillium and save into json
    if not (genus_term is None):
        Gen.ncbi_download(
            email=email,
            genus_term=genus_term,
            additional_term=additional_term,
            name_out=name_out,
            path_work=path_work,
            path_tmp=path_tmp,
            max_len=max_len,
        )
    elif not (accession_file is None):
        try:
            with open(accession_file, "r") as f:
                accession_list = [l.strip() for l in f.readlines()]
            Gen.ncbi_downloadbyacclist(
                email=email,
                list_acc=accession_list,
                name_out=name_out,
                path_work=path_work,
                path_tmp=path_tmp,
                max_len=max_len,
            )
        except:
            Gen.Mes(f"Accession file {accession_file} is not valid!")
            raise Exception
    else:
        Gen.Mes("Either genus name nor accession file is needed!")
        raise Exception

    Gen.Mes("GenMine finished. Please site us as below")
    Gen.Mes(
        "Chang Wan Seo, Sung Hyun Kim, Young Woon Lim & Myung Soo Park (2022) Re-Identification on Korean Penicillium Sequences in GenBank Collected by Software GenMine, Mycobiology, DOI: 10.1080/12298093.2022.2116816 https://www.tandfonline.com/doi/full/10.1080/12298093.2022.2116816"
    )
