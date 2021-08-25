# GenMiner
 A GenBank data mining program for (mostly fungal) taxonomists
 
# Install


# Use
python GenMine.py -e <your-email> -g <genus> -o <outpath>
 
Essential Parameters
--genus, -g : genus you want to parse
--out, -o : output path
 
Highly recommended Parameters
--email, -e : your email for NCBI access

Optional Parameters
--additional, -a : additional terms (ex. country name) to filter 
--max, -m : maximum length of the sequence to parse (default: 5000)
--start, -s : starting date of the record to parse
--end, -e : end date of the record to parse
 
# Features

 GenMiner is a python program that parses records from GenBank and sort by gene names, based on Entrez library.
 Comparing to Entrez, GenMiner has some advantages and disadvantages
 
 # Advantages
 - GenMiner doesn't misses records, especially with multiple terms
 - GenMiner can download discontinuously, especially useful in low internet condition
 - Classify downloaded records

# Limitations
- Slower than Entrez, due to completeness and stability
