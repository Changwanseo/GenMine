# GenMine_pypi/GenMine/command.py
import sys
import argparse
from importlib.metadata import version


class CommandParser:
    def __init__(self) -> None:
        # Argument input for running in commandline mode
        self.parser = argparse.ArgumentParser(
            description="Download genbank records", prog="GenMine"
        )

    def get_args(self) -> argparse.Namespace:
        self.parser.add_argument(
            "--email", "-e", help="Entrez Email, required", type=str, required=True
        )
        self.parser.add_argument(
            "--genus",
            "-g",
            nargs="*",
            help="List of genus to find | File with genera in each line",
            type=str,
        )
        self.parser.add_argument(
            "--accession",
            "-c",
            nargs="*",
            help="List of accessions to get | File with accessions in each line",
            type=str,
        )
        self.parser.add_argument(
            "--additional", "-a", nargs="*", help="additional terms", type=str
        )
        self.parser.add_argument(
            "--out",
            "-o",
            help="Out file path. Use previous result directory name for continue run",
            type=str,
        )
        self.parser.add_argument(
            "--max",
            "-m",
            help="Maximum length of sequence to search, in order to remove genomic sequences",
            type=int,
        )

        self.parser.add_argument(
            "--version",
            action="version",
            version=f"GenMine {version('GenMine')}",
        )

        """
        self.parser.add_argument(
            "--start",
            "-s",
            help="Starting date of finding records, as term of YYYY-MM-DD",
            type=str,
        )
        self.parser.add_argument(
            "--end",
            "-f",
            help="Final date of finding records, as term of YYYY-MM-DD",
            type=str,
        )
        """

        return self.parser.parse_args()
