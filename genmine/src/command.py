# GenMine_pypi/GenMine/command.py
import argparse


class CommandParser:
    def __init__(self) -> None:
        # Argument input for running in commandline mode
        self.parser = argparse.ArgumentParser(
            description="Download genbank records", prog="GenMine"
        )

    def get_args(self) -> argparse.Namespace:

        self.parser.add_argument(
            "--email", "-e", help="Entrez Email, required", type=str
        )
        self.parser.add_argument(
            "--taxa", "-t", "--genus", "-g",
            dest="taxa",
            nargs="*",
            help="List of taxa to find | File with taxa in each line",
            type=str
        )
        self.parser.add_argument(
            "--accession", "-c", nargs="*", help="List of accessions to get | File with accessions in each line", type=str
        )
        self.parser.add_argument(
            "--additional", "-a", nargs="*", help="additional terms", type=str
        )
        self.parser.add_argument("--out", "-o", help="Out file path. Use previous result directory name for continue run", type=str)
        self.parser.add_argument(
            "--max",
            "-m",
            help="Maximum length of sequence to search, in order to remove genomic sequences",
            type=int,
        )
        self.parser.add_argument(
            "--api-key",
            "-k",
            help="NCBI API key for higher rate limits (10 req/s instead of 3)",
            type=str,
        )
        self.parser.add_argument(
            "--check-network",
            action="store_true",
            default=False,
            help="Run network diagnostics for NCBI connectivity and exit",
        )
        
        '''
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
        '''


        return self.parser.parse_args()
