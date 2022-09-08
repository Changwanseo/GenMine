# GenMine
 A GenBank data mining program for (mostly fungal) taxonomists
 
 GenMine downloads GenBank nucleotide data.
 GenMine filters downloaded data with frequently used genes in taxonomy.

 
Citation: 
Chang Wan Seo, Sung Hyun Kim, Young Woon Lim & Myung Soo Park (2022) Re-Identification on Korean Penicillium Sequences in GenBank Collected by Software GenMine, Mycobiology, DOI: 10.1080/12298093.2022.2116816
 https://www.tandfonline.com/doi/full/10.1080/12298093.2022.2116816
 
 
# Install
WIP (we are going to provide pip and conda installation options)

* currently possible method

1. Download this repository
2. Enter to the folder and run following commands
```
conda env create --file GenMine.yaml
conda activate GenMine
```


# Use

* Usage
```
python GenMine.py -e <your-email> -g <genus> -o <outpath>
```

* Essential Parameters
```
--genus, -g : genus (or higher taxa) you want to parse
--out, -o : output path
```
 
* Highly recommended Parameters
```
--email, -e : your email for NCBI access
```
* Optional Parameters
```
--additional, -a : additional terms (ex. country name) to filter 
--max, -m : maximum length of the sequence to parse (default: 5000)
--start, -s : starting date of the record to parse
--end, -e : end date of the record to parse
```

# Features

 GenMine is a python program that parses records from GenBank and sort by gene names, based on Entrez library.
 Comparing to Entrez, GenMiner has some advantages and disadvantages
 
 # Advantages
 - GenMine doesn't misses records, especially with multiple terms
 - GenMine can download discontinuously, especially useful in low internet condition
 - GenMine classifies downloaded records by gentypes (ITS, LSU, SSU, *BenA* etc...)
 * If you want more gene types, issue it!

# Limitations
- Slower than Entrez (sometimes a lot), due to completeness and stability

# Bug reports and Suggestions
- Bug reports and suggestions are available in Issues or directly to wan101010@snu.ac.kr
- However, we want GenMine to remain as small tool. For suggestions that are little further to purpose of GenMine might be accepted in our upcomming softwares
