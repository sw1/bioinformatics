library(stringr)
library(readr)
library(ggplot2)


############################
############################
############ R #############
############################
############################

a <- 5
b = 5

a==b
!(a==b)

5==5 & 5==5
5==5 & 5==6
5==5 | 5==6

5==6 & notavariable
5==6 && notavariable


v <- vector(mode='double',length=5)
length(v)
v
v[1]
v[2] <- 5
v
v + 5

w <- 1:5
w

w+v
w*v

w %*% v

vector(mode='numeric',length=5)
vector(mode='character',length=5)
vector(mode='integer',length=5)
vector(mode='logical',length=5)

v <- vector(mode='logical',length=5)
v
!v
as.integer(v)
as.integer(!v)

A <- array(0,dim=c(5,5,5))
str(A)
length(A)
dim(A)

A <- matrix(0,nrow=5,ncol=5)
rownames(A) <- c('a','b','c','d','e')
colnames(A) <- c('a','b','c','d','e')
str(A)

L <- list(first=1,second=c(1,2,3,4,5),something=A,yeah='apple')
str(L)

L['new'] <- c('some','stuff')
str(L)

DF <- data.frame(idx=1:5,h=c(53,42,45,69,54),name=c('a','c','v','d','e'))
str(DF)

DF$somenewstuff <- c(4,2,3,4,5)

group <- c('dog','cat','crow','dog','dog','cat')
group <- as.factor(group)
group
str(group)

group[1] <- 'apple'

group <- c('dog','cat','crow','dog','dog','cat')
table(group)
table(group[1:2])

group <- as.factor(group)
table(group)
table(group[1:2])

plot(1:5,1:5)
plot(1:5,1:5,col='red')

qplot(1:5,1:5,color='red') + geom_line() + geom_abline(slope=1.25,color='blue')


if (1 == 1) print('hello')
if (1 != 5) print('hello')

a <- 5
if (1 == a){
  print('hello')
}else{
  print('goodbye')
}

if (1 == a){
  print('hello')
}else if (5 == a){
  print('hello again')
}else{
  print('goodbye')
}


for (i in 1:5){
  print(i)
}

animals <- c('cat','dog','crow')
for (animal in animals){
  print(animal)
}

for (i in 1:length(animals)){
  print(animals[i])
}

seq_along(animals)
for (i in seq_along(animals)){
  print(animals[i])
}

animals <- NULL

length(animals)
1:length(animals)
for (i in 1:length(animals)){
  print(animals[i])
}

for (i in seq_along(animals)){
  print(animals[i])
}


############################
############################
### FASTA AND BIOSTRINGS ###
############################
############################

library(reutils)

source('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/162b0c542482d481f79b0160071114eb38cb568e/r_bioinformatics_functions.R')


fasta <- readDNAStringSet('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/f1fb586160d12c34f29532c731066fd8912a0e0c/example.fasta',format='fasta')
fasta

View(fasta)

length(fasta)
width(fasta)
names(fasta)

subseq(fasta,start=1,end=5)
subseq(fasta,start=5,end=25)
subseq(fasta,start=5,width=10)
subseq(fasta,start=5,end=-5)

alphabetFrequency(fasta)
alphabetFrequency(fasta,baseOnly=TRUE)
alphabetFrequency(fasta,baseOnly=TRUE,as.prob=TRUE)
round(alphabetFrequency(fasta,baseOnly=TRUE,as.prob=TRUE),3)

letterFrequency(fasta,'G',as.prob=TRUE)
letterFrequency(fasta,'GC',as.prob=TRUE)
letterFrequency(fasta,'ACGT',as.prob=TRUE)
letterFrequency(fasta,'ACGT',OR=0,as.prob=TRUE)

consensusMatrix(fasta,baseOnly=TRUE)
consensusMatrix(fasta,baseOnly=TRUE,width=100)
dinucleotideFrequency(fasta)
trinucleotideFrequency(fasta)
oligonucleotideFrequency(fasta,width=5)

RNAStringSet(fasta)
reverseComplement(fasta)
translate(reverseComplement(fasta))




############################
############################
########## FASTQ ###########
############################
############################

fastq <- readDNAStringSet('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/f1fb586160d12c34f29532c731066fd8912a0e0c/example.fastq',format='fastq')
fastq




############################
############################
## CREATING SEQUENCE SETS ##
############################
############################

problem_createsequencesets()
s1
s2
s3

S <- c(s1,s2,s3)
SS <- DNAStringSet(S)
SS

names(SS) <- c('sequence_1','sequence_2','sequence_3')
SS

############################
############################
########### PASTE ##########
############################
############################


DOG <- c('dog1','dog2','dog3')
CAT <- c('cat1','cat2','cat3')

paste(DOG,CAT)
paste(DOG,CAT,sep='-')
paste(DOG,CAT,sep='_')
paste(DOG,CAT,sep='')
paste('dog','cat',1:3,sep='')
paste('dog','cat',1:3,sep='_')
paste('dog_','cat',1:3,sep='')

seq_names <- paste('sequence_',1:length(SS),' | User_12 | ',date(), sep='')
seq_names

names(SS) <- seq_names
SS
View(SS)

output_name <- 'seq_set_out.fasta'
writeXStringSet(SS,file=output_name,format="fasta")


############################
############################
######## METADATA ##########
############################
############################


FASTA <- readDNAStringSet('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/10bc2f50d1c739827ea2ba4edb146b36a6a4c14a/problems_metadata.fasta',format='fasta')
problem_metadata(FASTA)
FASTA
head(META)


which(META$Center == 'Philadelphia')
FASTA['Sequence_6333']

header_names <- META$ID[c(12,15,78)]
FASTA[header_names]

SS <- FASTA[META$Center == 'Houston' & META$Genus == 'Escherichia']
SS

output_name <- 'seq_subset_out.fasta'
writeXStringSet(SS,file=output_name,format="fasta")


############################
############################
####### GC FUNCTIONS #######
############################
############################


gc_calc <- function(x) (x['g']+x['c'])/sum(x)

gc <- function(s,pos){
  
  s <- stringr::str_to_lower(s)
  s <- unlist(strsplit(s,''))
  
  if (!missing(pos)) s <- s[seq(pos,length(s),3)]
  counts <- table(s)
  
  gc_calc(counts)
  
}


gc_skew_calc <- function(x) {counts <- table(x); (counts['g']-counts['c'])/(counts['g']+counts['c'])}

gc_skew <- function(s,win){
  
  s <- stringr::str_to_lower(s)
  s <- unlist(strsplit(s,''))
  
  if (missing(win)) {
    gc <- gc_skew_calc(s)
  }else{
    start <- seq(1,length(s),win)
    gc <- NULL
    for (i in start){
      gc <- c(gc, gc_skew_calc(s[(i):(i+win-1)]))
    }
  }
  
  gc
  
}


generate_random_dna_gc_s(len=1000,seed=1)

gc(s)
gc(s,1)
gc(s,2)
gc(s,3)

generate_random_dna_skew_s(len=1000,w=1,seed=5)

gc_skew(s)
gc_skew(s,100)
plot_skew(gc_skew(s,25))


############################
############################
######## NCBI SEARCH #######
############################
############################

ids1 <- esearch("CFTR AND human[Organism] AND complete",db='nucleotide',retmax=15,sort='relevance')
ids2 <- esearch("PKD1 AND human[Organism] AND complete",db='nucleotide',retmax=15,sort='relevance')
ids3 <- esearch("DMPK AND human[Organism] AND complete",db='nucleotide',retmax=15,sort='relevance')

ids_df <- reutils::content(esummary(ids1),'parsed')
View(ids_df)

efetch(ids1[1], rettype = "fasta", retmode = "text")
efetch(ids2[4], rettype = "fasta", retmode = "text")
efetch(ids3[5], rettype = "fasta", retmode = "text")

ids <- c(ids1[1],ids2[4],ids3[5])
ids

FASTA <- efetch(ids,db='nucleotide', rettype = "fasta", retmode = "xml")
LENS <- FASTA$xmlValue('//TSeq_length')
LENS
SEQS <- FASTA$xmlValue('//TSeq_sequence')
SEQS

tmp <- tempfile()
FASTA <- efetch(ids,db='nucleotide', rettype = "fasta", retmode = "text", outfile=tmp)
FASTA <- readDNAStringSet(tmp)






###################
#### SRA FILES ####
###################

library(SRAdb)

sqlfile <- '~/SRAdb/SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)

rs <- listSRAfile(c('SRP040765'), sra_con, fileType = 'sra')
str(rs)

run <- rs$run[1]
run

getSRAfile(run, sra_con, fileType = 'sra',destDir='~')

run_info <- getSRA(search_terms='SRP040765', out_types=c('run'),sra_con)
str(run_info)

sample_info <- getSRA(search_terms='SRP040765', out_types=c('sample'),sra_con)
str(sample_info)
sample_info$sample_attribute[1]

experiment_info <- getSRA(search_terms='SRP040765', out_types=c('experiment'),sra_con)
str(experiment_info)




###################
###### BLAST ######
###################


library(Biostrings)

S <- c('TGAAAAAGGCGAGCTGGTGGTTCTGGGACGCAACGGTTCCGACTACTCCGCT
        GCGGTGCTGGCGGCCTGTTTACGCGCCGATTGTTGCGAGATCTGGACGGATGTTGACGGTGTTTATACCT
        GCGATCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTATCAGGAAGCGATGGAGCTTTC
        TTACTTCGGCGCTAAAGTTCTTCACCCCCGCACCATCACCCCCATCGCCCAGTTTCAGATCCCTTGCCTG
        ATTAAAAATACCGGAAATCCTCAAGCACCAGGTACGCTCATTGGTGCCAGCCGTGATGAAGACGAATTAC
        CGGTCAAGGGCATTTCCAATCTGAATAACATGGCAATGTTCAGCGTTTCCGGCCCGGGGATGAAAGGGAT
        GGTTGGCATGGCGGCGCGCGTCTTTGCAGCGATGTCACGCGCCCGTATTTCCGTGGTGCTGATTACGCAA
        TCATCTTCCGAATACAGTATCAGTTTCTGCGTTCCGCAAAGCGACTGTGTGCGAGCTGAACGGGCAATGC
        AGGAAGAGTTCTACCTGGAACTGAAAGAAGGTTTACTGGAGCCGTTGGCGGTGACGGAACGGCTGGCCAT',
       'GTCAGAACCACGGGAAAATATCGTTTATCAGTGCTGGGAGC
        GTTTTTGCCAGGAGCTTGGCAAGCAAATTCCAGTGGCGATGACTCTGGAAAAGAATATGCCGATCGGTTC
        GGGCTTAGGCTCCAGCGCCTGTTCGGTGGTCGCGGCGCTGATGGCGATGAATGAACACTGCGGCAAGCCG
        CTTAATGACACTCGTTTGCTGGCTTTGATGGGCGAGCTGGAAGGACGTATCTCCGGCAGCATTCATTACG
        ACAACGTGGCACCGTGTTTTCTTGGTGGTATGCAGTTGATGATCGAAGAAAACGACATCATCAGCCAGCA
        AGTGCCAGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGGAATTAAAGTCTCGACGGCAGAAGCC
        CGGGCTATTTTACCGGCGCAGTATCGCCGCCAGGATTGCATTGCGCACGGGCGACATCTGGCTGGCTTCA
        TTCACGCCTGCTATTCCCGTCAGCCTGAGCTTGCCGCGAAGCTGATGAAAGATGTTATCGCTGAACCCTA
        CCGTGAACGGTTACTGCCTGGCTTCCGGCAGGCGCGGCAGGCGGTCGCGGAAATCGGCGCGGTAGCGAGC
        GGTATCTCCGGCTCCGGCCCGACCTTGTTCGCTCTATGTGACAAGCCGGATACCGCCCAGCGCGTTGCCG
        ACTGGTTGGGTAAGAACTACCTGCAAAATCAGGAAGGTT')
S <- gsub('[[:space:]]|\n','',S)
S <- DNAStringSet(S)
names(S) <- c('s1','s2')

FASTA <- tempfile()
writeXStringSet(S,file=FASTA,format="fasta")

SEQS <- readr::read_lines(FASTA)
SEQS

blastn <- '/data/sw1/ncbi-blast-2.5.0+/bin/blastn'

output <- c('sseqid', 'evalue', 'bitscore')

test <- system('pwd')
test

test <- system2('pwd',stdout=TRUE)
test

BLAST <-  system2(blastn,
                  args=c('-db',"'nt'",
                         '-outfmt',sprintf("'6 %s'",paste(collapse=' ',output)),
                         '-perc_identity',"'.99'",
                         '-entrez_query',"'Escherichia[ORGANISM]'",
                         '-remote'),
                  input=SEQS,
                  stdout=TRUE)

BLAST_OUT <- read.table(textConnection(BLAST),quote='',sep='\t',col.names=output)
head(BLAST_OUT)





###################
#### PHYLOSEQ #####
###################

library(phyloseq)
library(ape)
library(DESeq2)

OTU <- read.csv('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/1b05f24f189f14ea9902ac3867aca40c80ac6db3/otu_table.csv')
TAX <- read.csv('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/1b05f24f189f14ea9902ac3867aca40c80ac6db3/tax_table.csv')
SAMP <- read.csv('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/052dfdc3df97589f6405d79889c9b3b651eb1967/sample_metadata.csv')
TREE <- read.tree('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/052dfdc3df97589f6405d79889c9b3b651eb1967/tree.tree')

all(colnames(OTU) == SAMP$Sample_ID)
rownames(SAMP) <- SAMP$Sample_ID

TAX <- as.matrix(TAX)
rownames(TAX) <- paste0('otu',1:nrow(TAX))
rownames(OTU) <- rownames(TAX)

taxa_names(TREE) <- rownames(TAX)

PS <- phyloseq(otu_table(OTU,taxa_are_rows=TRUE),tax_table(TAX),sample_data(SAMP),phy_tree(TREE))

PS1 <- prune_samples(!is.na(sample_data(PS)$Enterotype),PS)
sample_data(PS1)$ENTEROTYPE <- as.factor(sample_data(PS1)$Enterotype)

PS1 <- filter_taxa(PS1,function(x) sum(x) > 0,prune = TRUE)

PS1

plot_tree(PS1,color='ENTEROTYPE')
plot_bar(PS1,fill='Group')