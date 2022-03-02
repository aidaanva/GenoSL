# GenoSL
Script to genotype single stranded libraries based on forward and reverse read mapping


## Usage 
```bash
Rscript genoSL.R [--] [--help] [--opts OPTS] [--output OUTPUT] [--fasta
       FASTA] input table
```
## Author:
Aida Andrades Valtueña (aida.andrades[at]gmail.com

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

Note: The names specified for each samples must correspond to the column name in the snpTable.tsv produced by MultiVCFAnalyzer
