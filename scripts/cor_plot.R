library(GGally)
require(circlize)
data<-read.csv("..//Data/my_data_match.csv",header=TRUE)

CM<-cor(data[,4:12])

ggpairs(data[,4:12])

ggcorr(data[,4:12], nbreaks = 10,label = TRUE)

require(psych)

op <- par(mfrow=c(3,2))
spider(y=1,x=2:9,data=Thurstone,connect=FALSE) #a radar plot
spider(y=1,x=2:9,data=Thurstone) #same plot as a spider plot
spider(y=1:3,x=4:9,data=Thurstone,overlay=TRUE)
#make a somewhat oversized plot
spider(y=1,x=2:9,data=CM) 
par(op)


chordDiagram(CM,symmetric = TRUE)


devtools::install_github("mattflor/chorddiag")
library(chorddiag)

## example taken from the github site
m <- matrix(c(11975,  5871, 8916, 2868,
              1951, 10048, 2060, 6171,
              8010, 16145, 8090, 8045,
              1013,   990,  940, 6907),
            byrow = TRUE,
            nrow = 4, ncol = 4)
haircolors <- c("black", "blonde", "brown", "red")
dimnames(m) <- list(have = haircolors,
                    prefer = haircolors)
m
#             prefer
#   have     black blonde brown  red
#     black  11975   5871  8916 2868
#     blonde  1951  10048  2060 6171
#     brown   8010  16145  8090 8045
#     red     1013    990   940 6907

groupColors <- c("#000000", "#FFDD89", "#957244", "#F26223")
chorddiag(CM)
