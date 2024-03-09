### Analysis of Clonal Expansion. 

### This script uses the Changeo-package shown in reference below: 
#   Gupta NT*, Vander Heiden JA*, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH. Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv359

library(devtools)
library(shazam)
library(alakazam)

### Step 1 ------------------------------------------------------------------
### For each patient, convert the heavy chain (H), lambda chain (L) and kappa chain (K) files (currently in .txz format and .fasta formats) into .tsv format for further analysis with Change-O. 
# List of all patients in this study
patient_list = c() #Contents of list have been deleted here to remove any patient identifiers from this code file. 

#patients 102048, 101106, 101164, 101005, 101104, 101220 and 102054, 102037, 101120 are excluded becuase the former two had no sequences isolated and the latter three had only 1 which makes analysis difficult

for (patient in patient_list) {
  name = as.character(patient) 
  H_chain_file_txz <- paste(name, "_H.txz", sep = "")
  H_chain_file_fasta <- paste(name, "_H.fasta", sep = "")
  H_chain_command <- paste("MakeDb.py imgt -i ", H_chain_file_txz, " -s ", H_chain_file_fasta, " --extended", sep = "")
  L_chain_file_txz <- paste(name, "_L.txz", sep = "")
  L_chain_file_fasta<- paste(name, "_L.fasta", sep = "")
  L_chain_command <- paste("MakeDb.py imgt -i ", L_chain_file_txz, " -s ", L_chain_file_fasta, " --extended", sep = "")
  K_chain_file_txz <- paste(name, "_K.txz", sep = "")
  K_chain_file_fasta <- paste(name, "_K.fasta", sep = "")
  K_chain_command <- paste("MakeDb.py imgt -i ", K_chain_file_txz, " -s ", K_chain_file_fasta, " --extended", sep = "")
  system(H_chain_command)
  system(L_chain_command)
  system(K_chain_command)
}

# In addition, create a file that will contain all the H_chain sequences from all patients combined. "merged.tsv" file. 
filenames_tsv_all <- as.character("")
for (patient in patient_list) {
  name = as.character(patient) 
  H_chain_file_tsv <- paste(name, "_H_db-pass.tsv", sep = "")
  filenames_tsv_all <- paste(filenames_tsv_all, " ", H_chain_file_tsv, sep = "")
}
filenames_tsv_all
command <- paste("ParseDb.py merge -d", filenames_tsv_all, " -o merged.tsv", sep = "")
system(command)


### Step 2: -----------------------------------------------------------------
### Run TigGer to obtain correct genotype calls. 
library(tigger)
library(dplyr)
library(readr)

# Load the germline databases that are newly downloaded - the SampleGermlineIgHV is a 2014 database. 
setwd("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/Human_IMGT_database")
germline_IGHV <- readIgFasta("IGHV.fasta")

# Load the merged file with all patients and find novel alleles
setwd("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/IMGT/Zip_and_fasta")
merged <- read_tsv("merged.tsv")
novel <- findNovelAlleles(merged,germline_IGHV, nproc=1, germline_min = 10, auto_mutrange = TRUE) 
selectNovel(novel)
plotNovel(merged, selectNovel(novel))

# Infer genotypes
geno <- inferGenotype(merged, germline_db=germline_IGHV, novel=novel,
                      find_unmutated=FALSE)
genotype_db <- genotypeFasta(geno, germline_IGHV, novel)
print(geno)
plotGenotype(geno, text_size = 10)

# Reassign the alleles. Updated genotype will be placed in the v_call_genotyped column
merged <- reassignAlleles(merged, genotype_db)

not_in_genotype <- merged$v_call %>%
  strsplit(",") %>%
  unlist() %>%
  unique() %>%
  setdiff(names(genotype_db))

data.frame(Ambiguous=c(mean(grepl(",", merged$v_call)),
                       mean(grepl(",", merged$v_call_genotyped))),
           NotInGenotype=c(mean(merged$v_call %in% not_in_genotype),
                           mean(merged$v_call_genotyped %in% not_in_genotype)),
           row.names=c("Before", "After")) %>% 
  t() %>% round(3)

write_csv(merged, "~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/IMGT/Zip_and_fasta/Post_tigger/merged_post_tigger.csv")


# All of the above code was for the merged H chain file. To assign correct genotypes to individual patient files, the code below subsets the merged dataframe into per patient dataframe 
for (patient in patient_list) {
  name = as.character(patient) 
  H_chain_file_tsv <- paste(name, "_H_db-pass.tsv", sep = "")
  file <- read_tsv(H_chain_file_tsv) 
  list_of_sequences <- file$sequence_id
  subset_merged <- merged[merged$sequence_id %in% list_of_sequences, ]
  directory <- paste("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/IMGT/Zip_and_fasta/Post_tigger/", H_chain_file_tsv, sep = "")
  write.table(subset_merged, directory, row.names = FALSE, sep = "\t")
}


### Step 3------------------------------------------------------------------
### Identify cutoff through Shazam.
### Step 1 has generated tsv files for all patients and chains. It has also generated .tsv file called merged.tsv that contains the combined H chains from patients with more than 1 H chain sequence isolated (>1, not >=1)

setwd("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/IMGT/Zip_and_fasta/Post_tigger")
library(readxl)
library(shazam)
library(readr)

# The below code determines a hamming distance cutoff value for all patients combined. 
# Set vCallColumn = v_call_genotyped becasue we have to let it know we inferred genotype
# Set first = FALSE for more leniency for ambiguous V genes

merged <- read_excel("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/IMGT/Zip_and_fasta/Post_tigger/merged_post_tigger.xlsx")
dist_ham_all_patients <- distToNearest(merged, sequenceColumn="junction", 
                                       vCallColumn="v_call_genotyped", jCallColumn="j_call",
                                       model="ham", normalize="len", nproc=1, first = FALSE)

library(ggplot2)
p1 <- ggplot(subset(dist_ham_all_patients, !is.na(dist_nearest)),
             aes(x=dist_nearest)) +
  theme_bw() +
  xlab("Hamming distance") +
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.2, color="firebrick", linetype=2)
plot(p1)

output <- findThreshold(dist_ham_all_patients$dist_nearest, method="density")
threshold <- output@threshold

# Plot distance histogram, density estimate and optimum threshold
plot(output, title="Density Method") 

### Step 4------------------------------------------------------------------
# Clonal explansion analysis using DefineClones.py. 
# The output of this section will be filename_clone-pass.tsv
setwd("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/IMGT/Zip_and_fasta/Post_tigger")

# Run DefineClones.py for each of the H chain files (files containing all H chains). Distance set to 0.2 as determined from the above method. 
system("DefineClones.py -d 1001_H_db-pass.tsv --act set --model ham    --norm len --dist 0.2 --vf v_call_genotyped")


for (patient in patient_list) {
  name = as.character(patient) 
  H_chain_file <- paste(name, "_H_db-pass.tsv", sep = "")
  H_chain_command <- paste("DefineClones.py -d ", H_chain_file, " --act set --model ham    --norm len --dist 0.2 --vf v_call_genotyped", sep = "")
  system(H_chain_command)
}


### Step 5 ------------------------------------------------------------------
### Count clones per patient 

# The output of DefineClones.py script does not contain the timepoint and subject id columns. The below code adds that information to the output of DefineClones.py.
library(alakazam)
library(writexl)
setwd("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq")
all_patients <- read_excel('BCR_seq_Medicago.xlsx') #This is a file that contains all of the sequence ids and their corresponding patient ID, timepoint and other information. 

# create new dataframe with all sequence ids and timepoints
timepoints <- all_patients[,c(1, 3)] #select columns 1 and 3 as those represent sequence_id and timepoint respectively. 


# Create a function that can take the last index of a string. This will be helpful for the section of code below it. 
lastIndexOf <- function(string) {
  last_index = nchar(string)
  character = substr(string, last_index, last_index)
  return(character)
}

# Create a function that matches and merges two dataframes by their Sequence IDs, data2 should be the timepoint dataframe
addTimepoint <- function(data1, data2) {
  data1 <- as.data.frame(data1)
  data2 <- as.data.frame(data2) 
  overlapping_ids = intersect(data1$`sequence_id`, data2$`sequence_id`)
  new_timepoints <- data2[data2$`sequence_id` %in% overlapping_ids, ] #subset the timepoint dataframe 
  result <- merge(data1, new_timepoints, by = "sequence_id", all = TRUE)
  return(result)
}

# Running the for loop for all patient files. This for loop will take the clone-pass.tsv files and add to them their proper timepoints.
# it will then output an excel file for all patients and store them in the folder Clones_by_patient folder
setwd("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/IMGT/Zip_and_fasta/Post_tigger")

for (patient in patient_list) {
  name = as.character(patient) 
  file_directory <- as.character("~/Documents/Pulendran_Lab/Medicago/Data/BCRseq/Patients_Fasta/IMGT/Zip_and_fasta/Post_tigger/Clones_by_patient/")
  filename_tsv = paste(name, "_H_db-pass_clone-pass.tsv", sep="")
  patient_clones <- read_tsv(filename_tsv)
  patient_clones_timepoint <- addTimepoint(patient_clones, timepoints)
  patient_summary <- countClones(patient_clones_timepoint, groups = "Timepoint", clone = "clone_id")
  outputfilename_xlsx = paste(file_directory, name, ".xlsx", sep="")
  write_xlsx(patient_summary, outputfilename_xlsx)
}

# The output of this section is a file containing clone counts for each patient. 

