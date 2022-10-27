# GenMine_pypi/GenMine/command.py
import argparse


class CommandParser:
    def __init__(self) -> None:
        # Argument input for running in commandline mode
        self.parser = argparse.ArgumentParser(
            description="Download genbank records", prog="GenMine"
        )

    def get_args(self) -> argparse.Namespace:

        self.parser.add_argument("--email", "-e", help="entrez email", type=str)
        self.parser.add_argument(
            "--genus", "-g", nargs="*", help="genus to find", type=str
        )
        self.parser.add_argument(
            "--additional", "-a", nargs="*", help="additional terms", type=str
        )
        self.parser.add_argument(
            "--accession", "-c", help="Accession list file to get", type=str
        )
        self.parser.add_argument("--out", "-o", help="out path", type=str)
        self.parser.add_argument(
            "--max",
            "-m",
            help="maximum length of sequence to search, in order to remove genomic sequences",
            type=int,
        )
        self.parser.add_argument(
            "--start",
            "-s",
            help="starting date of finding records, as term of YYYY-MM-DD",
            type=str,
        )
        self.parser.add_argument(
            "--end",
            "-f",
            help="final date of finding records, as term of YYYY-MM-DD",
            type=str,
        )

        return self.parser.parse_args()
