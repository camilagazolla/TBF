> ###   Q1
<br>

This representation is used to inspect the read quality scores. The median score in each position (red line) is shown, as well as the quartiles of the score distribution (yellow bars), and the mean score (blue line). The quality decays especially at the end, which is common when sequencing using an Illumina platform. 

This a very good sequence data once the median quality scores per base is above 30 in all positions. 
If it was a genome sequence data, for example, I would straightly proceed to assembly. If it was the forward reads of a paired-end metabarcoding experiment, I would trim the first 10 nucleotides and truncate the reads at position 270 (trimming the last 30 nucleotides), being sure that this proceed would preserve the overlap between the reads.

I generally use the plotQualityProfile function of the DADA2 R package insted of FastQC because it work better for automations using R. 
Here is an example of the use of plotQualityProfile below. I have added some [example FASTQ files](https://github.com/camilagazolla/Trial-By-Fire/tree/main/example_data) here for reproducibility. 😀

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
<img src="https://github.com/camilagazolla/Trial-By-Fire/blob/main/F3D143_S209_L001_R1_001_quality.png" width="50%" height="50%">
Figure. Quality profiles using plotQualityProfile function.

<br><br>

>###   Q2
<br>

To complete the task I have employed fuctions from the Biostrings R package. This package contains classes and functions for representing biological strings such as DNA, RNA and amino acids. The result table is available at the link bellow.

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

# sort by Length
df <-df[order(df$Lenght),]

# write the tab delimited table
write.table(df, "read_counts.tab", sep = "\t", row.names = FALSE, quote = FALSE)
``` 
<br><br>

>###   Q3
<br>

FASTA and FASTQ are both text archives that hold sequence data and metadata. The FASTQ files are generally used for storing the data from NGS experiments, while the FASTA files are generally used to store reference data.
 
Differently from FASTA, FASTQ has a standardized way of encoding quality scores for each nucleotide, which is important to create the plot displayed on Question 1, for example. In addition, the quality scores are also important for the DADA2 pipeline that creates a parametric error model for every amplicon dataset in order to differentiate sequencing errors from real biological variation.

<br><br>

>###   Q4
<br>

A total of 34 different 10 base DNA strings were obtained from the total 1733 sequences. The most frequent was TGGTGTCAGC	(59.66%, 1034 counts), followed by TGGTGCCAGC (664 counts), TGTGTCAGCC (3 counts) and TGTGCCAGCC (2 counts).

<br>

📌 Task

- Extract the first 10 bases of a sample and count the strings

Can consume 🐡

- [CVM382 FASTA file](https://github.com/camilagazolla/Trial-By-Fire/blob/main/CVM382.fasta)

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

<br><br>

>###   Q5
<br>

I have aligned a complete 16S rRNA sequence (J01859.1) with some reads using EMBOSS merger ([result here!](https://github.com/camilagazolla/Trial-By-Fire/blob/main/EMBOSS_output.txt)). The amplicon represents the 16S rRNA V4 region because i) it presents a length of 300 bp, and ii) it is aligned next to the position of the widely used primers 515F and 806R. The 16S rRNA V4 region is a popular target for metabarcoding of Bacteria, however, Archaea sequences could also be present.

<br><br>


>###   Q6
<br>

A “diversity index” is a mathematical measure of the diversity of species in a community.  The  “alpha diversity” represents the diversity applicable to a single sample, comprising “richness” and “evenness ”.  

- Richness or “how many microorganisms?” : the number of different species
- Evenness or “are the microorganisms balanced to each other?”: describes if there is uniformity of species (equal abundance)

I have been employing only alpha diversity indices that are towards evenness measures (i.e. Simpson e, Pielou e, Camargo, Bulla's index, Evar index, Heip e and Strong, based on [Hagerty et al., 2020](https://doi.org/10.1371/journal.pone.0229204)) when analyzing microbiome datasets, because it provides more information about the community composition than richness. For example, consider two communities of 100 individuals each and composed of 10 different species. One community has 10 individuals of each species, the other has one individual of each of the nine species and 91 individuals of the tenth species. Which community is more diverse/uniform? Clearly the first, but both communities would have the same species richness. Another problem is the presence of singletons (a read that has been detected exactly once). There is a difficulty in differentiating a “rare error” from a real singleton, which exists in nature. Therefore, high levels of singletons make richness estimates wrong.

<br><br>


>###   Q7
<br>

- Chao1 is a a measure of richness, which gives more weight to rare species. As previouly explained, the presence of sequencing errors make this measure less suitable for characterizing communities.

- Shannon is referred to a metric that capture both richness and relative evenness. However, this metric was described as loaded strongly onto the richness according to [Hagerty et al., 2020](https://doi.org/10.1371/journal.pone.0229204).


<br><br>


>###   Q8
<br>


In order to test the correlation between the two columns I confirmed that the variables presented a linear relationship. Afterwards, the normality test and normal probability plot showed that the data are not normally distributed and should be correlated using non-parametric methods. The Spearman and Kendall statistical tests were perfomed in a loop that created a scatter plot and a text file with the test results.

<img src="https://github.com/camilagazolla/Trial-By-Fire/blob/main/spearman_corrplot.png" width="40%" height="40%">
Figure. Scatter plot displaying the Spearman test results.
<br><br>

According to the results obtained, there is a statistically significant (p-value < 2.2e-16) strong positive relationship between (rho = 0.95, tau = 0.84) *Lactobacillus* and Axis.1 values.  

Importantly, in this question the data was treated as if it was any number, since there are strong dependencies between the relative abundances of different ASVs/OTUs, characterizing the compositional structure of the data ([Gloor et al., 2017](https://doi.org/10.3389/fmicb.2017.02224)), making this test inappropriate. 

<br>

📌 Task

- Perform correlation between relative abundance data of two variables

Can consume 🐡

- [Table with relative abundances](https://github.com/camilagazolla/Trial-By-Fire/blob/main/test.csv)

Can provide 🍣

- [Scatter plot](https://github.com/camilagazolla/Trial-By-Fire/blob/main/spearman_corrplot.png)
- [.txt with statistical test details](https://github.com/camilagazolla/Trial-By-Fire/blob/main/spearman_results.txt)

<br>

**On R:**

``` 
library(ggplot2)
library(ggpubr)

# path to folder containing .csv file
path <- "C:/Users/camila/Desktop/Trial_by_fire/"

setwd(path)

# read the .csv
test <- read.csv("test.csv")

# scatter plot to see if the covariation is linear
ggplot(test, aes(x=Axis.1, y=Lactobacillus)) + geom_point() # =>  yes! linear!

# shapiro-Wilk normality test and normal probability plot
shapiro.test(test$Axis.1) 
ggplot(test, aes(sample=Axis.1))+stat_qq()

shapiro.test(test$Lactobacillus) 
ggplot(test, aes(sample=Lactobacillus))+stat_qq()

# => data are not normally distributed, using non-parametric correlation

# create scatter plot and test using Spearman and Kendall
for (i in c("spearman", "kendall")) {
  
  # scatter plot
  p <- ggscatter(test, x = "Axis.1", y = "Lactobacillus",
                 add = "reg.line",  # regressin line
                 add.params = list(color = "blue", fill = "lightgray"), # customize reg. line
                 conf.int = TRUE # add confidence interval
  ) + font("ylab", face = "italic") + # italics on genus name
    stat_cor(method = i)
  
  png(paste0(i,"_corrplot.png")) # save plot as .png
  plot(p) 
  dev.off()
  
  sink(paste0(i,"_results.txt")) # save the test results in a .txt
  print(cor.test(test$Axis.1,test$Lactobacillus, method=i))
  sink()
  closeAllConnections()
}

``` 

<br><br>








