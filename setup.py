# hey-insutance/setup.py
from setuptools import setup, find_packages

setup(
    name="GenMine",
    version="1.0.6.0",
    description="GenBank data miner for fungal taxonomists",
    author="Changwan Seo",
    author_email="wan101010@snu.ac.kr",
    url="https://github.com/Changwanseo/GenMine",
    python_requires=">= 3.8",
    packages=find_packages(include=["GenMine"]),
    install_requires=["biopython", "pandas", "openpyxl", "xlrd", "numpy", "xmltodict"],
    zip_safe=False,
    # important part
    entry_points={"console_scripts": ["GenMine = GenMine.main:main"]},
    package_data={},
    include_package_data=True,
    license="GNU GENERAL PUBLIC",
)
