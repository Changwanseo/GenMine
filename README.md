# GenMine
 A GenBank data mining program for (mostly fungal) taxonomists
 
 GenMine downloads GenBank nucleotide records
 GenMine filters downloaded data with frequently used genes in taxonomy.

 
Citation: 
Chang Wan Seo, Sung Hyun Kim, Young Woon Lim & Myung Soo Park (2022) Re-Identification on Korean Penicillium Sequences in GenBank Collected by Software GenMine, Mycobiology, DOI: 10.1080/12298093.2022.2116816

https://www.tandfonline.com/doi/full/10.1080/12298093.2022.2116816
  
 
## Install
* pip
```
pip install GenMine
```

* conda 
```
conda install -c cwseo GenMine
```

* If problem occurs please uninstall then install with specifying versions
```
pip uninstall GenMine
pip install GenMine==1.0.6
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

* Download data accession numbers
```
GenMine -e wan101010@snu.ac.kr -c ON417149.1 ON417150.1
```

### Advanced usage
* Download records of multiple genera
```
GenMine -e wan101010@snu.ac.kr -g Penicillium Trichoderma Alternaria
```

* Download records of multiple genera given by file
```
GenMine -e wan101010@snu.ac.kr -g genera.txt
```
"genera.txt" should be like this
```
Penicillium
Trichoderma
Alternaria
```

* Download records of multiple accession given by file
```
GenMine -e wan101010@snu.ac.kr -c accessions.txt
```
"accessions.txt" should be like this
```
ON417149.1
ON417150.1
MW554209.1
OK643788.1
```

* Continue download from interrupted run (only for accessions, for genus, it will automatically solve if you launch GenMine in same location)
```
GenMine -e wan101010@snu.ac.kr -c accessions.txt -o "2022-11-02-00-12-08"
# Caution 1: -o should be name of previous run result directory
# Caution 2: will not work for finished run
```


* Parameters
```
--genus, -g : List of genus to find | File with genera in each line
--accession, -c : List of accessions to get | File with accessions in each line
--email, -e : your email for NCBI access
```
* Optional Parameters
```
--additional, -a : additional terms (ex. country name) to filter 
--max, -m : maximum length of the sequence to parse (default: 5000)
```

## Output explanations
### Main output

WIP

## Features

 GenMine is a python program that parses records from GenBank and sort by gene names, based on Entrez library.
 Comparing to Entrez, GenMiner has some advantages and disadvantages
 
 ### Advantages
 - GenMine doesn't misses records, especially with multiple terms
 - GenMine can download discontinuously, especially useful in low internet condition
 - GenMine classifies downloaded records by gene types (ITS, LSU, SSU, *BenA* etc...)
 
 * If you want more gene types, issue it!
 * We are currently working on better gene annotations

### Limitations
- Slower than Entrez (sometimes a lot), due to completeness and stability

## Bug reports and Suggestions
- Bug reports and suggestions are available in Github Issues or directly to wan101010@snu.ac.kr
- However, we want GenMine to remain as small tool. For suggestions little bit too much for the purpose of GenMine might be accepted in our upcomming softwares
