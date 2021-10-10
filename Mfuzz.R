library(Mfuzz)
library(MSnbase)
library(MSqRob)
library(biobroom)

pepdat <- import2MSnSet('peptides.txt', filetype="MaxQuant")




pepdat <- filterZero(pepdat, 0.25) # Filter zeroes

pepdat.r <- filter.NA(pepdat, 0.25)
pepdat.f <- fill.NA(pepdat.r, mode='mean')
tmp <- filter.std(pepdat.f, min.std=0)
pepdat.s <- standardise(pepdat.f)

cl <- mfuzz(pepdat.s, c=10, m=1.25)
mfuzz.plot(pepdat.s, cl=cl, mfrow=c(4,4))

m1 <- mestimate(pepdat.s)

cl2 <- mfuzz(pepdat.s, c=18, m=1.35)
mfuzz.plot(pepdat.s, cl=cl2, mfrow=c(4,4))

Dmin(pepdat.s,m1,crange=seq(2,40,1),repeats=3,visu=TRUE)

acore.list <- acore(pepdat.s, cl=cl2)
temp   = sapply(acore.list,'[[',"NAME")

pepdat.s@featureData@data$Leading.razor.protein
