library(readr)
library(Biostrings)
library(rentrez)

############################
############################
########## BLAST ###########
############################
############################

########## REMOTE ##########

S <- c('ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA',
       'TCAGACCTGAGTGGCGCTAACCATCCGGCGCAGGCAGGCGATTTGCAGTACGGCTGGAATTGTCAC')
S <- DNAStringSet(S)
names(S) <- c('s1','s2')

FASTA <- tempfile()
writeXStringSet(S,file=FASTA,format="fasta")

SEQS <- readr::read_lines(FASTA)
SEQS



blastn <- '/data/sw1/ncbi-blast-2.5.0+/bin/blastn'


test <- system('pwd')
test

test <- system2('pwd',stdout=TRUE)
test

test <- system2('pwd',stdout=FALSE)
test

# https://www.ncbi.nlm.nih.gov/books/NBK279675/

output <- c('sseqid', 'evalue', 'bitscore')

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

BLAST_OUT <- tempfile()
write_lines(BLAST,path=BLAST_OUT)
read_delim(BLAST_OUT,delim='\t',col_names=output)



########## LOCAL ##########

ids <- entrez_search("Escherichia[ORGANISM]",db='nucleotide',retmax=2500,
                     use_history=TRUE) # required if retmax>500

n <- 100
FASTA <- tempfile()
for (i in seq(1,ids$retmax,n)){
  seqs <- entrez_fetch(db='nucleotide', rettype = "fasta", retmode = "text", 
                       web_history=ids$web_history,retstart=i,retmax=n)
  cat(seqs,file=FASTA,append=TRUE)
  cat('Appended ',i+n-1,' of ',ids$retmax,' sequences.\n',sep='')
}
DB_SEQS <- readr::read_lines(FASTA)
head(DB_SEQS)

makeblastdb <- '/data/sw1/ncbi-blast-2.5.0+/bin/makeblastdb'
db_name <- '/data/sw1/Bioinformatics/db_ecoli/db_ecoli'

system2(makeblastdb,
        args=c('-out',db_name,
               '-dbtype','nucl',
               '-title','ecoli',
               '-parse_seqids',
               '-hash_index'),
        input=DB_SEQS,
        stdout=FALSE)

output <- c('sseqid', 'evalue', 'bitscore')

BLAST <-  system2(blastn,
                  args=c('-db',db_name,
                         '-outfmt',sprintf("'6 %s'",paste(collapse=' ',output)),
                         '-perc_identity',"'.99'"),
                  input=SEQS,
                  stdout=TRUE)

BLAST_OUT <- read.table(textConnection(BLAST),quote='',sep='\t',col.names=output)
head(BLAST_OUT)


############################
############################
######### ggplot2 ##########
############################
############################

library(dplyr)
library(ggplot2)
library(babynames)
library(ggrepel)

set.seed(123)

N <- 25
DAT <- data.frame(height=c(rnorm(N,68,3),rnorm(N,64,2)),
                  weight=c(rnorm(N,170,8),rnorm(N,125,5)),
                  sex=rep(c('M','F'),each=N),
                  name=sample(babynames$name,size=N*2,replace=TRUE))
                  
ggplot() + geom_point(data=DAT,aes(x=height,y=weight))

ggplot() + geom_point(data=DAT,aes(x=height,y=weight,color=sex))

ggplot() + geom_point(data=DAT,aes(x=height,y=weight,color=sex)) + 
  geom_hline(yintercept=c(125,170)) +
  geom_vline(xintercept=c(64,68)) 

ggplot() + geom_point(data=DAT,aes(x=height,y=weight,color=sex)) + 
  geom_rug(data=DAT,aes(x=height,y=weight,color = sex))

ggplot() + geom_point(data=DAT,aes(x=height,y=weight),color='black') + 
  geom_rug(data=DAT,aes(x=height,y=weight,color=sex)) +
  stat_density2d(data=DAT,aes(x=height,y=weight,fill=sex,alpha=..level..), geom = "polygon")

ggplot() + 
  geom_rug(data=DAT,aes(x=height,y=weight,color=sex)) +
  stat_density2d(data=DAT,aes(x=height,y=weight,fill=sex,alpha=..level..), geom = "polygon") +
  geom_point(data=DAT,aes(x=height,y=weight),color='black')

ggplot() + 
  geom_rug(data=DAT,aes(x=height,y=weight,color=sex)) +
  stat_density2d(data=DAT,aes(x=height,y=weight,fill=sex,alpha=..level..), geom = "polygon") +
  geom_point(data=DAT,aes(x=height,y=weight),color='black') +
  guides(alpha=FALSE)

FIG <- ggplot() + 
  geom_rug(data=DAT,aes(x=height,y=weight,color=sex)) +
  stat_density2d(data=DAT,aes(x=height,y=weight,fill=sex,alpha=..level..), geom = "polygon") +
  geom_point(data=DAT,aes(x=height,y=weight),color='black') +
  guides(color=FALSE,alpha=FALSE)

FIG + ggtitle('Body Type')

last_plot() + labs(x='Height (inches)',y='Weight (lbs)',fill='Sex')

last_plot() + labs(x='Height (inches)',y='Weight (lbs)',fill='Sex') +
  geom_label_repel(data=DAT %>% filter(height>70 & weight>170),
                   aes(x=height,y=weight,label=name,fill=sex),
                   fontface='bold',color='black')


############################
############################
######### PHYLOSEQ #########
############################
############################

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

OTU <- otu_table(OTU,taxa_are_rows=TRUE)
TAX <- tax_table(TAX)
SAMP <- sample_data(SAMP)
TREE <- phy_tree(TREE)

PS <- phyloseq(OTU,TAX,SAMP,TREE)

PS1 <- prune_samples(!is.na(sample_data(PS)$Enterotype),PS)
sample_data(PS1)$ENTEROTYPE <- as.factor(sample_data(PS1)$Enterotype)

PS1 <- filter_taxa(PS1,function(x) sum(x) > 0,prune=TRUE)

PS
PS1

plot_richness(PS1,x='ENTEROTYPE',color='ENTEROTYPE')

plot_bar(PS1,fill='Group')

PS1RA <- transform_sample_counts(PS1, function(x) x/sum(x))
plot_bar(PS1RA,fill='Group')

ORD <- ordinate(PS1,method='MDS',distance='bray')
plot_ordination(PS1,ORD,color='ENTEROTYPE') + geom_point(size=5)
plot_ordination(PS1,ORD,color='ENTEROTYPE') + geom_point(size=5) + facet_wrap(~SeqTech)

ORD <- ordinate(PS1,method='NMDS')
plot_ordination(PS1,ORD,color='ENTEROTYPE') + geom_point(size=5) + facet_wrap(~SeqTech)

ORD <- ordinate(PS1,method='CCA')
plot_ordination(PS1,ORD,color='ENTEROTYPE') + geom_point(size=5) + facet_wrap(~SeqTech)

PS2 <- subset_taxa(PS1, Group=='Group1')
PS2 <- prune_samples(sample_data(PS2)$SeqTech == 'Pyro454',PS2)
PS2 <- filter_taxa(PS2,function(x) sum(x) > 0,prune=TRUE)
plot_tree(PS2,color='ENTEROTYPE',base.spacing=0.005,size='abundance')
plot_tree(PS2,color='ENTEROTYPE',base.spacing=0.005,size='abundance') + coord_polar(theta='y')

PS3 <- prune_samples(sample_data(PS1)$SeqTech == 'Pyro454',PS1)
taxa_filtered <- names(sort(taxa_sums(PS3),decreasing=TRUE))[1:50]
PS3 <- prune_taxa(taxa_filtered,PS3)
plot_heatmap(PS3,sample.label='ENTEROTYPE',sample.order='ENTEROTYPE',method='MDS',distance='bray',max.label=500)

NET <- make_network(PS3,max.dist=.3,distance='bray')
plot_network(NET,PS3,color='ENTEROTYPE',label=NULL)

PS4 <- prune_samples(!is.na(sample_data(PS)$Gender),PS)
PS4 <- filter_taxa(PS4,function(x) sum(x) > 0,prune=TRUE)
diagdds <- phyloseq_to_deseq2(PS4, ~ Gender)
diagdds <- DESeq(diagdds, test='Wald', fitType='parametric')
res <- results(diagdds, cooksCutoff=FALSE)
res


BIOM <- import_biom('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/93eac154ed375a3e37fe7a76fdf2419e98280222/biom_file.biom')
?import_biom
?parse_taxonomy_default

# https://joey711.github.io/phyloseq/import-data.html#_microbio_me_qiime