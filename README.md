# GenoSL
Script to genotype single stranded libraries based on forward and reverse read mapping

# Installation
Requirements: have a working installation of conda. 

Clone this github directory. 

Create conda environment by running:
```bash
conda env create -f environment.yaml
```
Activate conda environment:
```bash
conda activate genoSL
```
Run the script as shown in usage.

## Usage 
```bash
Rscript genoSL.R [--] [--help] [--opts OPTS] [--output OUTPUT] [--fasta
       FASTA] input table
```
## Author:
Aida Andrades Valtueña (aida.andrades[at]gmail.com )

## Description:
```bash
Genotypes single stranded libraries based on forward and reverse
mapping reads to account for potential damage

positional arguments:
  input         Path to the snpTable.tsv produced by MultiVCFAnalyzer
  table         Path to table.tsv containing genomes to Genotype. See
                README.md for specific formating

flags:
  -h, --help    show this help message and exit

optional arguments:
  -x, --opts    RDS file containing argument values
  -o, --output  Specify root of output name, default is snpTable
                [default: snpTable]
```

## Input
Table: table.tsv containing genomes to Genotype.
The table should be formated as:
| sampleName | All              |	Reverse           | Forward            |
|------------|------------------|--------------------|--------------------|
| Sample1    |Sample1AllReadsVCF| Sample1ReverseReads| Sample1ForwardReads|
| Sample2    |Sample2AllReadsVCF| Sample2ReverseReads| Sample2ForwardReads|
| Sample3    |Sample3AllReadsVCF| Sample3ReverseReads| Sample3ForwardReads|

Note: The names specified for each samples must correspond to the column name in the snpTable.tsv produced by MultiVCFAnalyzer.

## Concept

In the absence of repair mechanism, we will observe an accumulation of DNA damage at the end of the fragments over time, which is one of the characteristics used to authenticate ancient DNA. One of the typical DNA damages observed is the deamination of cytosines (C) leading to uracils (U), which occurs due to hydrolytic damage. When constructing DNA libraries with fragments containing deaminated C, these are interpreted as T by the polymerase. In the case of single-stranded libraries (sslib) the damage appears only as C to T substitutions. This is due to the absence of blunt-end repairing during the library construction, which will result in the reads containing misscoding substitutions in the form of A to G substitutions in the 3’ end of reads. 

This characteristic of sslib can be used to determine if an observed change in a nucleotide in comparison with the reference is due to damage or not. Reads mapping in the forward direction (5’ to 3’) will present damage as C to T  substitutions, while reads mapping in the reverse direction (3’ to 5’) will present damage as G to A substitutions (this is due to the reverse complementing of the original read). In that regard, the strategy employed here will be that any T call could be due to damage, so we will look at the reverse mapping reads to confirm if this is a true call or a false positive due to damage. With the same logic, A calls could be due to damage and we will look at the forward mapping reads in order to determine if it is a true call or a false call due to damage. 

Disclaimer: This method essentially halves the callable data for any T or A call.
