##############################################################################################################
#########   input - path to the FASTA file of DNA sequences of the same length
#########           regionLength
#########           binSize
#########           sortRegionSize
#########   output - jpg file with heatmap of GC content centered around middle of the provided seuqences
#########          - jpg file with plot of average GC content of provided sequences aligned at the center
##############################################################################################################

library(stringr) 
library(hexbin)
library(RColorBrewer)
library(stringi)
library(grid)
library(matrixStats)

#############  Variables - user input
file_name<-".\\worm_TES.fasta"
regionLength<-2000
binSize<-50
sortRegionSize<-0
numberOfHexbins<-100
##################################
#{
my_fasta<-readLines(file_name)
sss<-my_fasta[seq(2,length(my_fasta),2)]

############## Verify user input
if(regionLength != round(regionLength)){stop("Region length should be integer.")}
if(binSize != round(binSize)){stop("Bin size should be integer.")}
if(binSize >= regionLength){stop("Bin size should be less than region length.")}
if(length(unique(nchar(sss))) != 1){stop("Input sequences must be of the same length.")}
if(nchar(sss[1])<regionLength){stop("Region length must be equal or smaller than the length of the input sequences.")}
if(length(which(str_count(sss,"A|a|T|t|C|c|G|g") != nchar(sss))) != 0){stop("Sequences contain characters other than A, T, C and G.")}

############## replace A and T with 0, C and G with 1
sss_bin <- gsub("(A|T)", "0", sss, ignore.case = TRUE)
sss_bin <- gsub("(C|G)", "1", sss_bin, ignore.case = TRUE)

center_of_region<-round(nchar(sss[1])/2)
############## coerce strings to a dataframe
df<-as.data.frame(strsplit(sss_bin,split=""))
ii<-which(df[(center_of_region-round(regionLength/2)):(center_of_region+round(regionLength/2)),]==1, arr.ind=TRUE)

xx<-ii[,1]
yy<-ii[,2]
############### Normalize GC % so that the highest displayed velue is 60 %
############### and the lowest is 18 %. Hexbins with values outside of this range
############### will be displayed without shading
######################################################################################
max_count<-(ncol(df)*regionLength/(numberOfHexbins^2))*0.6
min_count<-(ncol(df)*regionLength/(numberOfHexbins^2))*0.18

################ Plot GC % histogram
bin<-hexbin(xx, yy, xbins=numberOfHexbins)
my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))

jpeg(".\\GC_Heatmap.jpg")
plot(bin,colramp=my_colors ,maxcnt=max_count,mincnt=min_count,ylab="",xlab="",main="GC content") 

#### Change text in axis and legend of the hexbin plot
grid.edit(grep("xaxis",grid.ls(print=FALSE)$name,value=TRUE),gp=gpar(fontsize=14),at = c(0,regionLength/2,regionLength),label=c(paste(-round(regionLength/2),"bp", sep=" "),"center",paste(round(regionLength/2),"bp", sep=" ")))
grid.edit(grep("yaxis",grid.ls(print=FALSE)$name,value=TRUE),gp=gpar(fontsize=14))
my_label<-c("","18","20","23","25","28","30","33","36","39","42","45","48","50","52","55","58","60","gc %")
grid_text_elem<-grep("text",grid.ls(print=FALSE)$name,value=TRUE)
for(i in 2:19){grid.edit(grid_text_elem[i],label=my_label[i])
}
dev.off()

###Plot average GC profile
y<-rowSums(sapply(df[(center_of_region-round(regionLength/2)):(center_of_region+round(regionLength/2)),],as.numeric))/ncol(df)
x<-seq(0,length(y)-1,1)

bx<-seq(0,regionLength,binSize)
y_binned<-binMeans(y, x = x, bx = bx)

jpeg(".\\GC_AveragePlot.jpg")
xx<-seq(1,length(y_binned),1)
plot(xx,y_binned*100,type="l",main="Average GC %", ylab="% GC",xlab="",xaxt="n",cex.axis=1.5,cex.lab=1.5,lwd=3)
abline(v=(1+(length(xx)-1)/2),lty=2)
axis(1, at=c(xx[1],(1+(length(xx)-1)/2),max(xx)),cex.axis=1.5, labels=c(paste(-round(regionLength/2),"bp", sep=" "),"center",paste(round(regionLength/2),"bp", sep=" ")))
dev.off()
#}

