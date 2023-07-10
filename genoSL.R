#!/usr/bin/env Rscript
library(argparser)
library(tidyverse)
library(seqinr)
library(data.table)


parser <- argparser::arg_parser( 'Genotypes single stranded libraries based on forward and reverse mapping reads to account for potential damage',
                                 name = 'genoSL.R')
parser <- add_argument(parser, 'input',
                       type='character',
                       nargs=1,
                       help='Path to the snpTable.tsv produced by MultiVCFAnalyzer')
parser <- add_argument(parser, 'table',
                       type = 'character',
                       nargs = 1,
                       help = 'Path to table.tsv containing genomes to Genotype. See README.md for specific formating')
parser <- add_argument(parser, '--output',
                       type = 'character',
                       help = 'Specify root of output name, default is snpTable',
                       default = 'snpTable')
parser <- add_argument(parser, '--fasta',
                       type = 'character',
                       help = 'Indicate a fullAlignment.fasta to include the corrected genotyped bases for the sample',
                       default = 'None')
argv <- parse_args(parser)


##FUNCTIONS:

'%!in%' <- function(x,y)!('%in%'(x,y))

genotyping <- function(df, forward, reverse, aReads, sampleName){
  snpTable <- df 
  names<-colnames(snpTable)
  names <- replace(names, names==forward, "forwardReads")
  names <- replace(names, names==reverse, "reverseReads")
  names <- replace(names, names==aReads, "allReads")
  colnames(snpTable) <- names
  
  snpTableGenotyped <- snpTable %>%
    mutate(genotypedSample=case_when(
      reverseReads == "T" | forwardReads == "T" | allReads == "T" ~ reverseReads,
      forwardReads == "A" | reverseReads == "A" | allReads == "A" ~ forwardReads,
      TRUE ~ allReads))
  
  correctGenotype <- paste(sampleName, "c", sep = "_")
  names<-colnames(snpTableGenotyped)
  names <- replace(names, names=="forwardReads", forward)
  names <- replace(names, names=="reverseReads", reverse)
  names <- replace(names, names=="allReads", aReads)
  names <- replace(names, names=="genotypedSample", correctGenotype)
  colnames(snpTableGenotyped) <- names
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
writingFasta <- function(df, output){
  df2 <- df %>%
    gather(Samples,Call, 3:ncol(.)) %>%
    mutate(CallT = ifelse(Call==".", Ref, Call)) %>%
    select(-Ref, -Call) %>%
    spread(Position, CallT)
    
  write_tsv(df2, output, col_names = F)
  outputFasta <- paste(output, "fasta", sep = ".")
  system(paste("awk '{print \">\"$1; $1=\"\"; print $0}' OFS=", output, ">", outputFasta))
  lsystem(paste("rm",output))
  print("Fasta_dp saved")
}
fastatoTibble <- function(fasta) {
  fastaList <- seqinr::read.fasta(fasta, whole.header = T, forceDNAtolower = F)
  fastaM <- data.table(matrix(unlist(fastaList),
                             nrow = length(fastaList),
                             byrow = T))
  fastaM$Genome <- names(fastaList)
  df_sampleRows <- fastaM %>%
    pivot_longer(-contains("Genome"), names_to = "Position",values_to = "Call") %>%
    mutate(Position = gsub("V", "", Position)) #%>%
    #select(-PositionIn)
  return(df_sampleRows)
  print("fastaToTibble Comple")
}

### END FUNCTIONS

print("Reading snpTable and genoTable")
snpTable <- read.delim(argv$input, sep = "\t", stringsAsFactors = F, check.names = F)
genoTable <- read.delim(argv$table, sep = "\t", stringsAsFactors = F, check.names = F) %>% 
  mutate(sampleCorrected = paste(sampleName, "c", sep = "_"))
snpTableGT <- snpTable

print ("snpTable and genoTable done. Genotyping started")
for (row in 1:nrow(genoTable)) {
  nameForwardReads <- genoTable[row, "Forward"]
  nameReverseReads <- genoTable[row, "Reverse"]
  allReads <- genoTable[row, "All"]
  sampleName <- genoTable[row, "sampleName"]
  snpTableGT <- genotyping(snpTableGT, nameForwardReads, nameReverseReads, allReads, sampleName)
#  snpTable_c <- snpTableGT %>%
#    select(-nameForwardReads, -nameReverseReads, -allReads) %>%
#    preparingAllTable(sampleName)
  
#  snpTable_AllIncluded <- preparingAllTable(snpTableGT, sampleName)
}

print("Genotyping finished")
#To do if a reference fasta is provided:
allIncluded <- paste(argv$output, "allColumns.tsv", sep = "_")
corrected <- paste(argv$output, "genotyped.tsv", sep = "_")
fastaFile <- paste(argv$output, "genotyped", sep = "_")

snpTable_c <- snpTableGT[,!names(snpTableGT) %in% genoTable$All]
snpTable_c <- snpTable_c[,!names(snpTable_c) %in% genoTable$Reverse]
snpTable_c <- snpTable_c[,!names(snpTable_c) %in% genoTable$Forward]


if(argv$fasta != "None" ){
  print("fullAlignment.fasta has been provided, start run with fasta mode")
  #Extract sequence to genotype (need to Figure out)
  fullFasta <- seqinr::read.fasta(argv$fasta, whole.header = T, forceDNAtolower = F)
  print("fullFasta read")
  #Remove Forward and Reverse entries in the fasta
  fullFastamod <- fullFasta[names(fullFasta) %in% c(genoTable$Reverse,genoTable$Forward) == FALSE]
  #Select only genotyped data
  genotypedOnly <- snpTable_c %>%
    select(Position, Ref, genoTable$sampleCorrected) %>%
    mutate_at(vars(-Position, -Ref),
              funs(case_when(
                . == "." ~ Ref,
                TRUE ~ .)))
  #For loop to modify vectors within the fullFasta list
  "Include genotyping in fullFasta"
  for(i in genoTable$sampleName){
    name <- paste(i, "c", sep = "_")
    names(fullFastamod)[names(fullFastamod) == i] <- "name"
    posToChange <- genotypedOnly %>% select(Position, name)
    names <-colnames(posToChange)
    names <- replace(names, names==name, "name")
    colnames(posToChange) <- names
    fullFastamod$name[posToChange$Position] <- posToChange$name
    names(fullFastamod)[names(fullFastamod) == "name"] <- name
  }
  
  fullfastaFile <- paste(argv$output, "genotyped", "fullAlignment.fasta", sep = "_")
  seqinr::write.fasta(sequences = fullFastamod, names = names(fullFastamod), nbchar = 80, file.out = fullfastaFile)
  
  print("fullAlignment.fasta genotyped saved")
  
}

print("writing genotyped snpTable & snpFasta")

writingFasta(snpTable_c, fastaFile)

write_tsv(snpTableGT, allIncluded)
write_tsv(snpTable_c, corrected)

print("All done!")
