##VCF Preparation

library("vcfR")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Filtering of invariant sites from vcf

#####-------------------------------------------------------#####

##functions needed for filtering

replace_DNA<-function(matrix){
  matrix[matrix=="0/0"]<-"R"
  matrix[matrix=="1/1"]<-"A"
  matrix[matrix=="./0"]<-"R"
  matrix[matrix=="./1"]<-"A"
  matrix[matrix=="0/1"]<-"H"
  matrix[matrix==".:."]<-"N"
  matrix[matrix=="./."]<-"N"
  return(matrix)
}

variants<-function(df){
  invariant<-vector()
  r<-unique(df)
  r<-sort(r)
  if(length(r)==1){invariant<-c(invariant, TRUE)}
  else if(length(r)==2){
    if(all(r==c("A","R"))){invariant<-c(invariant, FALSE)}else{invariant<-c(invariant, TRUE)}
  }
  else if(length(r)==3){
    if(all(r==c("A","H","R")) | all(r==c("A","N","R"))){invariant<-c(invariant, FALSE)}else{invariant<-c(invariant, TRUE)} 
  }
  else if(length(r)==4){invariant<-c(invariant, FALSE)}
  return(invariant)
}

## load DATASET
vcf<-read.vcfR(file = "PI1.vcf")

## check DATASET
head(vcf)
head(vcf@gt)
head(vcf@fix)
head(vcf@fix[,1])                     ## get chromosom name (class = character)
head(vcf@fix[,2])                     ## get position of snp (class = character)

## remove invariant sites
df<-vcf@gt
df<-df[,-1]
df[]<-sapply(df, FUN = function(x) substr(x, 1, 3))

df2<-replace_DNA(df)
invariant<-apply(df2, MARGIN = 1, FUN = variants)

remove_inv<-which(invariant==T)
vcf@fix<-vcf@fix[-remove_inv,]
vcf@gt<-vcf@gt[-remove_inv,]

## SAVE new vcf file
write.vcf(vcf, file = "PI2.vcf.gz")

#####------------------------------------------------------------#####

## Additional filter options: minimal distance + subsample

#####------------------------------------------------------------#####
## filter by distance between loci
pos<-as.numeric(vcf@fix[,2])
dist<-c(0,diff(pos))
dist[which(dist<0)]<-0
remove_dist<-which(dist!=0 & dist<200)

vcf@fix<-vcf@fix[-remove_dist,]
vcf@gt<-vcf@gt[-remove_dist,]

## subsample DATASET
count<-c(1:length(vcf@fix[,2]))
set.seed(123)
selection<-sample(x = count, size = 50000)
selection<-sort(selection)

vcf@fix<-vcf@fix[selection,]
vcf@gt<-vcf@gt[selection,]

## SAVE new vcf file
write.vcf(vcf, file = "PI2_filtered.vcf.gz")

#####------------------------------------------------------------#####