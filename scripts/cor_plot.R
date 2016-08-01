library(GGally)
require(circlize)
require(gplots)
require(MASS)
data<-read.csv("..//Data/my_data_match.csv",header=TRUE)
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
CM<-cor(data[,4:12])

column_annotation <- sample(c("red", "blue", "green"), 9, replace=T)
column_annotation <- as.matrix(column_annotation)
colnames(column_annotation) <- c("Variable X")
row_annotation <- sample(c("red", "blue", "green"), 9, replace=T)
row_annotation <- as.matrix(t(row_annotation))
rownames(row_annotation) <- c("Variable Y")
heatmap.3(CM,RowSideColors=row_annotation, ColSideColors=column_annotation)
parcoord(data[,4:12])


# load packages and prepare data
library(alluvial)
tit <- as.data.frame(Titanic)

# only two variables: class and survival status
tit2d <- aggregate( Freq ~ Class + Survived, data=tit, sum)

alluvial( tit2d[,1:2], freq=tit2d$Freq, xw=0.0, alpha=0.8,
          gap.width=0.1, col= "steelblue", border="white",
          layer = tit2d$Survived != "Yes" )

heatmap.2(CM,trace="none",col=bluered(10))



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
