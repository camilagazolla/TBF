> ###   Q1
> **Interpret the plot below. What can be said about the overall quality of sequencing run used to generate this figure?**
> <img src="https://github.com/camilagazolla/Trial-By-Fire/blob/main/quality_scores_example.png" width="50%" height="50%">
<br>

This kind of representation is used to inspect the read quality scores. The median score in each position (red line) is shown, as well as the quartiles of the score distribution (yellow bars), and the mean score (blue line). The quality decays especially at the end, which is common when sequencing using an Illumina platform. I would consider this a very good sequence data. The trimmer parameters I would employ in this case would be trimming the first 10 nucleotides and the truncation of the reads at position ~270 (trimming the last 30 nucleotides).
Given the style of the plot, I believe it was done with the FastQC program. I used to employ this program a lot, however, today I generally use the plotQualityProfile function of the DADA2 R package, simply because it worked better for the automations I wanted to perform. 
I intend to show an example of the use of plotQualityProfile below. I have added some [example FASTQ files](https://github.com/camilagazolla/Trial-By-Fire/tree/main/example_data) here for reproducibility. 😀

<br>

📌 Task 
- Plot quality profiles for FASTQ files

Can consume 🐡
- [FASTQ files](https://github.com/camilagazolla/Trial-By-Fire/tree/main/example_data)

Can provide 🍣

- [Quality profile plots](https://github.com/camilagazolla/Trial-By-Fire/blob/main/F3D143_S209_L001_R2_001_quality.png)
<br>

**On R:**

``` 
# load package
library(dada2)

# path to folder containing the reads (decompressed FASTQ)
path <- "C:/Users/camila/Desktop/Trial_by_fire/example_data"

setwd(path)

# create "quality" folder on path
dir.create("quality") 

# file parsing
fns <- sort(list.files(pattern=".fastq"))

# plot quality profile for each FASTQ
for (i in unique(fns)){
  png(paste0(path, "/quality/", gsub(".fastq","", i), "_quality.png"))
  print(plotQualityProfile(i))
  dev.off()
}

``` 
<br><br>

>###   Q2
>**Generate a tab delimited table summarizing the read counts for all of the samples in the “fastas” folder. There should be a “SampleID” column and a read count column. Use what ever approach you prefer, but bonus point will be given for making a simple loop and using grep. Explain what you did and provide any accompanying code with comments.**

<br>

📌 Task

- Extract the number of sequences ("read conts") from a series of multi-FASTA files

Can consume 🐡

- FASTA files

Can provide 🍣

- [Table with read counts](https://github.com/camilagazolla/Trial-By-Fire/blob/main/read_counts.tab)
<br>

**On R:**

``` 
# load package for FASTA reading
library(Biostrings) 

# path to folder containing the FASTA files
path <- "C:/Users/camila/Desktop/Trial_by_fire/fastas/"

setwd(path)

# file parsing
fileNames <- sort(list.files(pattern=".fasta"))

# loop to extract "read counts" 
df  <-  NULL # crate df to store the results
for (i in fileNames){
  SampleID <- gsub(".fasta","",i) # remove extension from the name
  fas <- readDNAStringSet(i) # read the FASTA file
  Lenght <- length(fas) # extract how many sequences it contains
  df <-  rbind(df, data.frame(SampleID, Lenght)) # populating the df
}

# sort by Lenght
df <-df[order(df$Lenght),]

# write the tab delimited table
write.table(df, "read_counts.tab", sep = "\t", row.names = FALSE, quote = FALSE)
``` 
<br><br>

>###   Q3
>**What is a fasta? Is it different from fastq? If so, how?**
<br>

write later

<br><br>

>###   Q4
>**Get the first 10 bases of each read in sample “CVM382” found within the “fasta” directory. What 10 base strings repeats the most and what is its frequency? Explain what you did and provide any accompanying code with comments.**
<br>

A total of 34 different 10 base DNA strings were obtained from the total 1733 sequences. The most frequent was TGGTGTCAGC	(59.66%, 1034 counts), followed by TGGTGCCAGC (664 counts), TGTGTCAGCC (3 counts) and TGTGCCAGCC (2 counts).

<br>

📌 Task

- Extract the first 10 bases of a sample and count the strings

Can consume 🐡

- [CVM382 FASTA file] (https://github.com/camilagazolla/Trial-By-Fire/blob/main/CVM382.fasta)

Can provide 🍣

- [Table with string counts](https://github.com/camilagazolla/Trial-By-Fire/blob/main/CVM382Freq.tab)

<br>

**On R:**

``` 
library(Biostrings) 

# path to folder containing the FASTA files
path <- "C:/Users/camila/Desktop/Trial_by_fire/fastas/"

setwd(path)

# read the FASTA file
CVM382 <- readDNAStringSet("CVM382.fasta") 

# truncate at position 10
CVM382Subseq <- subseq(CVM382, end= 10) 

# create crosstab
CVM382Subseq.df <- as.data.frame(table(CVM382Subseq))

# sort by Freq
CVM382Subseq.df <- CVM382Subseq.df[order(CVM382Subseq.df$Freq),]

# write the tab delimited table
write.table(CVM382Subseq.df, "CVM382Freq.tab", sep = "\t", row.names = FALSE, quote = FALSE)
```








