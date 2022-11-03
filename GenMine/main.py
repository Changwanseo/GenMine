def main():

    from . import GenMine as Gen
    from GenMine.command import CommandParser
    from datetime import datetime, date

    import os

    args = CommandParser().get_args()

    # Managing arguments

    # default setup
    # date_start = -1
    # date_end = -1
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
        print("Emails are mandatory for safe NCBI connection")
        raise Exception

    # Get accession inputs
    accession_list = []
    if not (args.accession is None):
        if type(args.accession) is str:  # For only one input
            args.accession = [args.accession]
        for acc in args.accession:
            if os.path.isfile(acc):
                try:
                    with open(acc, "r") as f:
                        accession_list += [x.strip() for x in f.readlines()]
                except:
                    print(f"{acc} was recognized as input file, but cannot be parsed")
                    print(
                        f"The format of accession input file should be one accession in each line"
                    )
                    raise Exception
            else:
                accession_list.append(acc)

    # Filter accession_list
    accession_list = Gen.filter_acc(accession_list, email)

    # Get genus inputs
    genus_list = []
    if not (args.genus is None):
        if type(args.genus) is str:
            args.genus = [args.genus]
        for acc in args.genus:
            if os.path.isfile(acc):
                try:
                    with open(acc, "r") as f:
                        genus_list += [x.strip() for x in f.readlines()]
                except:
                    print(f"{acc} was recognized as input file, but cannot be parsed")
                    print(
                        f"The format of genus input file should be one genus in each line"
                    )
                    raise Exception
            else:
                genus_list.append(acc)

    # Abort when no inputs available
    if len(accession_list) == 0 and len(genus_list) == 0:
        print("No genus term nor accession numbers inserted. Aborted")
        raise Exception

    # Get additional terms, such as Korea
    if args.additional != None:
        additional_term = list(args.additional)

    # Get output location
    if args.out != None:
        path_out = args.out
    else:
        path_out = None

    # Get maximum sequence length
    if args.max != None:
        max_len = args.max

    # For time stamp
    date = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    # Download all seqs from NCBI with genus Penicillium and save into json
    if len(genus_list) > 0:
        for genus in genus_list:
            # Set output file name appendix
            name_out = genus

            # Set output location and file name
            if not (path_out is None):
                path_work = path_out
            else:
                path_work = f"{os.getcwd()}/{name_out}"

            # Set temporary directory location
            path_tmp = f"{path_work}/tmp"

            # Make running directory if not exists
            try:
                os.mkdir(path_work)
            except:
                pass

            # Make temporary directory if not exists
            try:
                os.mkdir(path_tmp)
            except:
                pass

            # Download data here
            Gen.ncbi_download(
                email=email,
                genus_term=genus,
                additional_term=additional_term,
                name_out=name_out,
                path_work=path_work,
                path_tmp=path_tmp,
                max_len=max_len,
            )

    if len(accession_list) > 0:
        # Set output file name appendix
        name_out = date

        # Set output location
        if not (path_out is None):
            path_work = path_out
        else:
            path_work = f"{os.getcwd()}/{name_out}"

        # Set temporary directory location
        path_tmp = f"{path_work}/tmp"

        # Make running directory if not exists
        try:
            os.mkdir(path_work)
        except:
            pass

        # Make temporary directory if not exists
        try:
            os.mkdir(path_tmp)
        except:
            pass

        # Download data here
        Gen.ncbi_downloadbyacclist(
            email=email,
            list_acc=accession_list,
            name_out=name_out,
            path_work=path_work,
            path_tmp=path_tmp,
            max_len=max_len,
        )

    Gen.Mes("GenMine finished. Please site us as below")
    Gen.Mes(
        "Chang Wan Seo, Sung Hyun Kim, Young Woon Lim & Myung Soo Park (2022) Re-Identification on Korean Penicillium Sequences in GenBank Collected by Software GenMine, Mycobiology, DOI: 10.1080/12298093.2022.2116816 https://www.tandfonline.com/doi/full/10.1080/12298093.2022.2116816"
    )
