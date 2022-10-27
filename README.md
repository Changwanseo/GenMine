# GenMine
 A GenBank data mining program for (mostly fungal) taxonomists
 
 GenMine downloads GenBank nucleotide data.
 GenMine filters downloaded data with frequently used genes in taxonomy.

 
Citation: 
Chang Wan Seo, Sung Hyun Kim, Young Woon Lim & Myung Soo Park (2022) Re-Identification on Korean Penicillium Sequences in GenBank Collected by Software GenMine, Mycobiology, DOI: 10.1080/12298093.2022.2116816
 https://www.tandfonline.com/doi/full/10.1080/12298093.2022.2116816
 
 
## Install
* pip
```
pip install GenMine
```

* conda (WIP)
Currently, please use indirect method

```
conda create -n GenMine
conda install pip
pip install GenMine
```

* If problem occurs please uninstall then install with specifying versions
```
pip uninstall GenMine
pip install GenMine==1.0.3.0
```

## Usage

### Basic usage

* Download all Penicillium records
```
GenMine -e wan101010@snu.ac.kr -g Penicillium
```

* Download all Korean penicillium records
```
GenMine -e wan101010@snu.ac.kr -g Penicillium -a Korea
```

* Download data with given accession numbers
```
GenMine -e wan101010@snu.ac.kr -c input.txt
```
"input.txt" should be like this

```
ON417149.1
ON417150.1
MW554209.1
OK643788.1
```


* Parameters
```
--genus, -g : genus (or higher taxa) you want to parse
--accession, -c : file with accession number in each row
--email, -e : your email for NCBI access
```
* Optional Parameters
```
--additional, -a : additional terms (ex. country name) to filter 
--max, -m : maximum length of the sequence to parse (default: 5000)
--start, -s : starting date of the record to parse
--end, -e : end date of the record to parse
```

## Features

 GenMine is a python program that parses records from GenBank and sort by gene names, based on Entrez library.
 Comparing to Entrez, GenMiner has some advantages and disadvantages
 
 ### Advantages
 - GenMine doesn't misses records, especially with multiple terms
 - GenMine can download discontinuously, especially useful in low internet condition
 - GenMine classifies downloaded records by gene types (ITS, LSU, SSU, *BenA* etc...)
 * If you want more gene types, issue it!

### Limitations
- Slower than Entrez (sometimes a lot), due to completeness and stability

## Bug reports and Suggestions
- Bug reports and suggestions are available in Issues or directly to wan101010@snu.ac.kr
- However, we want GenMine to remain as small tool. For suggestions that are little further to purpose of GenMine might be accepted in our upcomming softwares
