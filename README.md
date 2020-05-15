# GenoSL
Script to genotype single stranded libraries based on forward and reverse read mapping


## Usage 
```bash
Rscript genoSL.R [--] [--help] [--opts OPTS] [--output OUTPUT] [--fasta
       FASTA] input nameForwardReads nameReverseReads allReads
       sampleName
```
## Author:
Aida Andrades Valtueña (aida.andrades[at]gmail.com

## Description:
```bash
Genotypes single stranded libraries based on forward and reverse
mapping reads to account for potential damage

positional arguments:
  input             Path to the snpTable.tsv produced by
                    MultiVCFAnalyzer
  nameForwardReads  Name of the column in snpTable.tsv corresponding to
                    forward reads
  nameReverseReads  Name of the column in snpTable.tsv corresponding to
                    forward reads
  allReads          Name of the column in snpTable.tsv corresponding to
                    all reads
  sampleName        specify name of the sample, the new column in
                    snpTable and header in fasta file will be
                    sampleNamec

flags:
  -h, --help        show this help message and exit

optional arguments:
  -x, --opts        RDS file containing argument values
  -o, --output      Specify root of output name, default is snpTable
                    [default: snpTable]
  -f, --fasta       Include a fullAlignment.fasta to include the
                    corrected genotyped bases for the sample
```
