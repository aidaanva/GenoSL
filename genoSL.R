#!/usr/bin/env Rscript
library(argparser)
library(tidyverse)
library(seqinr)


parser <- argparser::arg_parser( 'Genotypes single stranded libraries based on forward and reverse mapping reads to account for potential damage',
                                 name = 'genoSL.R')
parser <- add_argument(parser, 'input',
                       type='character',
                       nargs=1,
                       help='Path to the snpTable.tsv produced by MultiVCFAnalyzer')
parser <- add_argument(parser, 'nameForwardReads',
                       type="character",
                       nargs=1,
                       help='Name of the column in snpTable.tsv corresponding to forward reads')
parser <- add_argument(parser, 'nameReverseReads',
                       type='character',
                       nargs=1,
                       help='Name of the column in snpTable.tsv corresponding to forward reads')
parser <- add_argument(parser, 'allReads',
                       type = 'character',
                       nargs=1,
                       help='Name of the column in snpTable.tsv corresponding to all reads')
parser <- add_argument(parser, 'sampleName',
                       type = 'character',
                       help = 'specify name of the sample, the new column in snpTable and header in fasta file will be sampleName_c')
parser <- add_argument(parser, '--output',
                       type = 'character',
                       help = 'Specify root of output name, default is snpTable',
                       default = 'snpTable')
parser <- add_argument(parser, '--fasta',
                       type = 'character',
                       help = 'Include a fullAlignment.fasta to include the corrected genotyped bases for the sample',
                       default = 'None')
argv <- parse_args(parser)


##FUNCTIONS:

'%!in%' <- function(x,y)!('%in%'(x,y))

genotyping <- function(df, forward, reverse, aReads, sampleName){
  snpTable <- read.delim(df, sep = "\t", stringsAsFactors = F)
  names<-colnames(snpTable)
  names <- replace(names, names==forward, "forwardReads")
  names <- replace(names, names==reverse, "reverseReads")
  names <- replace(names, names==aReads, "allReads")
  names <- str_replace(names, "^X", "")
  names <- str_replace(names, "^X", "XX")
  colnames(snpTable) <- names
  
  snpTableGenotyped <- snpTable %>%
    mutate(genotypedSample=case_when(
      reverseReads == "T" | forwardReads == "T" | allReads == "T" ~ reverseReads,
      forwardReads == "A" | reverseReads == "A" | allReads == "A" ~ forwardReads,
      TRUE ~ allReads))
    return(snpTableGenotyped)
}
preparingAllTable <- function(df, sampleName){
  snpTableGenotyped <- df
  names <- colnames(snpTableGenotyped)
  forwardReads <- paste(sampleName, "F", sep = "_")
  reverseReads <-paste(sampleName, "R", sep = "_")
  allReads <- paste(sampleName, "all", sep = "_")
  correctGenotype <- paste(sampleName, "c", sep = "_")
  
  names <- replace(names, names=="forwardReads", forwardReads)
  names <- replace(names, names=="reverseReads", reverseReads)
  names <- replace(names, names=="allReads", allReads)
  names <- replace(names, names=="genotypedSample", correctGenotype)
  colnames(snpTableGenotyped) <- names
  return(snpTableGenotyped)
}
fastatoTibble <- function(fasta) {
  fastaList <- seqinr::read.fasta(fasta, whole.header = T, forceDNAtolower = F)
  fastaM <- as_tibble(matrix(unlist(fastaList),
                             nrow = length(fastaList),
                             byrow = T),
                      .name_repair = "universal")
  fastaM$Genome <- names(fastaList)
  df_sampleRows <- fastaM %>%
    pivot_longer(-contains("Genome"), names_to = "PositionIn",values_to = "Call") %>%
    mutate(Position = gsub("...", "", PositionIn)) %>%
    select(-PositionIn)
  return(df_sampleRows)
  print("fastaToTibble Comple")
}

snpTableGT <- genotyping(argv$input, argv$nameForwardReads, argv$nameReverseReads, argv$allReads, argv$sampleName)
snpTable_c <- snpTableGT %>%
  select(-forwardReads, -reverseReads, -allReads) %>%
  preparingAllTable(argv$sampleName)

snpTable_AllIncluded <- preparingAllTable(snpTableGT, argv$sampleName)


if(argv$fasta != "None" ){
  print("fullAlignment.fasta has been provided, start run with fasta mode")
  #Extract sequence to genotype (need to Figure out)
  baseNameFullFasta <- paste(argv$output, sampleName, "fullAlignment.fasta", sep = "_")
  system(paste("awk '/^>/ {printf(\"\n%s\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\n\");}\' <'", args$fasta, "| grep -f id.txt -A 1 -m 1 >", baseNameFullFasta, sep = ""))
  fastaFull <- fastatoTibble("~/extracted.fasta")
}

allIncluded <- paste(argv$output, argv$sampleName, "allColumns.tsv", sep = "_")
corrected <- paste(argv$output, argv$sampleName, "c.tsv", sep = "")

write_tsv(snpTable_AllIncluded, allIncluded)
write_tsv(snpTable_c, corrected)

snpTableForFasta <-snpTable_c %>% gather(Samples,Call, 3:ncol(.))
snpTableForFasta$Call = ifelse(snpTableForFasta$Call==".", as.character(snpTableForFasta$Ref),as.character(snpTableForFasta$Call))
snpTable_c_fsnpalignment <- snpTableForFasta %>%
  spread(Samples, Call)
forFasta <- paste(argv$output, argv$sampleName,"c_forfasta.tsv", sep = "")
write_tsv(snpTable_c_fsnpalignment, forFasta)
