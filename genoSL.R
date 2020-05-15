#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<5) {
  stop("At least four argument must be supplied. Run script with:
       Rscript genotyping_script.r /path/to/multiVCFAnalyzer/results Name_Column_Forward_reads Name_Column_Reverse_reads Name_Column_all_reads Sample_name (output_Name)
       output_Name is optional, default is snpTable",
       call.=FALSE)
} else if (length(args)==5) {
  # default output file
  args[6] = "snpTable"
}

library(tidyverse)
##FUNCTIONS:

'%!in%' <- function(x,y)!('%in%'(x,y))


#input file
df <- paste(args[1])
#output file
baseNameOutput <- paste(args[6])

#Read snpTable
snpTable <- read.delim(df, sep = "\t", stringsAsFactors = F)
sampleName <- args[5]

print(args[2])
print(args[3])

names<-colnames(snpTable)
names <- replace(names, names==args[2], "forwardReads")
names <- replace(names, names==args[3], "reverseReads")
names <- replace(names, names==args[4], "allReads")
names <- str_replace(names, "^X", "")
names <- str_replace(names, "^X", "XX")
#print(names)
colnames(snpTable) <- names

genotypedSample <- paste(sampleName, "F", sep = "_")
snpTable_mod <- snpTable %>%
  mutate(genotypedSample=case_when(
    reverseReads == "T" | forwardReads == "T" | allReads == "T" ~ reverseReads,
    forwardReads == "A" | reverseReads == "A" | allReads == "A" ~ forwardReads,
    TRUE ~ allReads
  ))
#snpTable_mod %>% select(Ref, forwardReads, reverseReads, allReads, genotypedSample)

snpTable_c <- snpTable_mod %>%
  select(-forwardReads, -reverseReads, -allReads)

forwardReads <- paste(sampleName, "F", sep = "_")
reverseReads <-paste(sampleName, "R", sep = "_")
allReads <- paste(sampleName, "all", sep = "_")
correctGenotype <- paste(sampleName, "c", sep = "_")

names<-colnames(snpTable_mod)
names <- replace(names, names=="forwardReads", forwardReads)
names <- replace(names, names=="reverseReads", reverseReads)
names <- replace(names, names=="allReads", allReads)
names <- replace(names, names=="genotypedSample", correctGenotype)
colnames(snpTable_mod) <- names

names<-colnames(snpTable_c)
names <- replace(names, names=="forwardReads", forwardReads)
names <- replace(names, names=="reverseReads", reverseReads)
names <- replace(names, names=="allReads", allReads)
names <- replace(names, names=="genotypedSample", correctGenotype)
colnames(snpTable_c) <- names

allIncluded <- paste(baseNameOutput, sampleName, "allColumns.tsv", sep = "_")
corrected <- paste(baseNameOutput, sampleName, "c.tsv", sep = "")

write_tsv(snpTable_mod, allIncluded)
write_tsv(snpTable_c, corrected)

snpTableForFasta <-snpTable_c %>% gather(Samples,Call, 3:ncol(.))
snpTableForFasta$Call = ifelse(snpTableForFasta$Call==".", as.character(snpTableForFasta$Ref),as.character(snpTableForFasta$Call))
snpTable_c_fsnpalignment <- snpTableForFasta %>%
  spread(Samples, Call)
forFasta <- paste(baseNameOutput, sampleName,"c_forfasta.tsv", sep = "")
write_tsv(snpTable_c_fsnpalignment, forFasta)
