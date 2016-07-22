#library(devtools)
#install_github("ropensci/rebi", force=T)
#install_github("chrislad/phenotypicForest", force=T)
library(europepmc)
library(pubmed.mineR)
library(rebi)
library(stringr)
library(ggplot2)
library(plyr)
library(phenotypicForest)
library(sjPlot)
library(rcrossref)


abstracts <- readabs("pubmed_result.txt")
bio_oxford <- epmc_search(query='ISSN:1367-4803', limit = 40000) #Bio Oxf
bio_evol_online <- epmc_search(query='ISSN:1176-9343', limit = 40000) #Bio Evo Online
bio_bmc <- epmc_search(query='ISSN:1471-2105', limit = 40000) #BMC
gen_res <- epmc_search(query='ISSN:1088-9051', limit = 40000) #genome research
plos_compbio <- epmc_search(query='ISSN:1553-734X', limit = 40000) #Plos Comp Bio
jj <- rbind(bio_bmc, bio_oxford, bio_evol_online, plos_compbio, gen_res)


# let's start by searching for abstracts which mention common programming langugaes
langs <- searchabsL(abstracts, include=c("python","java","javascript",
                                         "perl","visual basic","cobol","c\\+\\+","c#",
                                         "php",'ruby',"matlab"))

# don't forget about 1 letter programming languages (C, R)
langs_one_let <- searchabsL(abstracts, include=c(" R ", " C "))

# now, let do some hacky parsing: for each abstract get PMID and programming languages
languages<-data.frame(matrix(nrow=length(langs@PMID),ncol=2))

# loop over all abstracts that we selected at a first place
for(i in 1:length(langs@PMID)){

  # let's keep track of where we are
  print(i)

  # store a PMID
  languages[i,1]<-(langs@PMID[i])

  # get rid of end of line symbols to have a single-line abstract for grep
  cc<-str_replace_all(langs@Abstract[i], "[\r\n]" , "")

  # now grep for keywords, ignore case
  cc<-str_extract_all(cc, regex("python|java|javascript|perl|visual basic|cobol|c\\+\\+|c#|php|ruby|matlab", ignore_case = TRUE))

  # sometimes more than 1 language is use, let's count which is a primary language for this abstract
  vv <- names(sort(table(cc[[1]]),decreasing=T))

  # if more than 1 - take the most frequently mentioned language
  if (length(vv) > 0){
    languages[i,2]<-tolower(vv[1])
  }
}

# now let's repeat the same for single letter programming languages
languages_one_let<-data.frame(matrix(nrow=length(langs_one_let@PMID),ncol=2))
for(i in 1:length(langs_one_let@PMID)){
  print(i)
  languages_one_let[i,1]<-(langs_one_let@PMID[i])
  cc<-str_replace_all(langs_one_let@Abstract[i], "[\r\n]" , "")
  cc<-str_extract_all(cc, regex(" R | C ", ignore_case = TRUE))
  vv <- names(sort(table(cc[[1]]),decreasing=T))
  if (length(vv) > 0){
    languages_one_let[i,2]<-tolower(vv[1])
  }
}

# bind both single letter and normal language matricies together
ll <- rbind(languages,languages_one_let)

# give meaningful names to columns
names(ll) <- c("pmid","lang")

# merge with journal data
dd <-merge(jj, ll, by="pmid",all.x = TRUE)


# finally let's plot the result with a polar Histogram to see the evolution of language popularity
kk <- (dd[!is.na(dd$lang),c(11,29)])
bb <- ddply(kk, .(pubYear), function(d) {data.frame(table(d$lang)/length(d$lang))})
colnames(bb) <- c("item","score","value")
bb$family <- "all"
pp <- polarHistogram(bb, familyLabel = FALSE, direction = "outwards")
print(pp +  scale_fill_manual(values= c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                         '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                                         '#cab2d6','#6a3d9a','#ffff99')))
