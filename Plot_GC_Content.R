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
library(optparse)

#############  Variables - user input
option_list = list(
make_option(c("-f", "--file"), type="character", default=NULL, help="sequences file name", metavar="character"),
                   make_option(c("-r", "--regionLength"), type="integer", default=1000, help="legnth of the region to be plotted", metavar="integer"),
                   make_option(c("-b", "--binSize"), type="integer", default=5, help="bin size to be used in average plot", metavar="integer"),
                   make_option(c("-l", "--label"), type="character", default="center", help="label to put on the plots", metavar="character"),
                   make_option(c("-s", "--sortRegionSize"), type="integer", default=0, help="size of the region to use for sorting; if 0 sequences remain unsorted", metavar="integer"),
                   make_option(c("-n", "--numberOfHexbins"), type="integer", default=100, help="number of hexbins to use in GC% heatmap; must be between 10 and 100", metavar="integer")
 );
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

out_file_heatmap<-paste(opt$file,"_heatmap.jpg",sep="")
out_file_plot<-paste(opt$file,"_plot.jpg",sep="")

regionLength<-opt$regionLength
binSize<-opt$binSize
label<-opt$label
sortRegionSize<-opt$sortRegionSize
numberOfHexbins<-opt$numberOfHexbins
##################################
{

my_fasta<- readLines(opt$file)
if(length(my_fasta) < 2) {stop("Bad fasta.")}
sss<-my_fasta[seq(2,length(my_fasta),2)]

############## Verify user input
if(numberOfHexbins != round(numberOfHexbins)){stop("Number Of Hexbins should be an integer.")}
if(numberOfHexbins < 10 || numberOfHexbins > 100){stop("Number Of Hexbins should be between 10 and 100.")}
if(sortRegionSize != round(sortRegionSize)){stop("Sort Region Size should be an integer.")}
if(sortRegionSize > round(regionLength)){stop("Sort Region Size should be less than region length.")}
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

##########################################################################
####### If sorting option is specified by the user,
####### sort sequences by the GC content of the specified window size
##########################################################################

if(sortRegionSize > 0){sort_value<-colSums(sapply(df[(center_of_region-round(sortRegionSize/2)):(center_of_region+round(sortRegionSize/2)),],as.numeric))
                       df_sorted<-df[order(sort_value)]
                       df<-df_sorted}


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

jpeg(out_file_heatmap,width=400, height=400)
plot(bin,colramp=my_colors ,maxcnt=max_count,mincnt=min_count,ylab="",xlab="",main="GC content") 

#### Change text in axis and legend of the hexbin plot
grid.edit(grep("xaxis",grid.ls(print=FALSE)$name,value=TRUE),gp=gpar(fontsize=14),at = c(0,regionLength/2,regionLength),label=c(paste(-round(regionLength/2),"bp", sep=" "),label,paste(round(regionLength/2),"bp", sep=" ")))
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

jpeg(out_file_plot,width=400, height=400)
xx<-seq(1,length(y_binned),1)
plot(xx,y_binned*100,type="l",main="Average GC %", ylab="% GC",xlab="",cex.axis=1.5,cex.lab=1.5,xaxt="n",lwd=2)
abline(v=(1+(length(xx)-1)/2),lty=2)
axis(1, at=c(xx[1],(1+(length(xx)-1)/2),max(xx)),cex.axis=1.5, labels=c(paste(-round(regionLength/2),"bp", sep=" "),label,paste(round(regionLength/2),"bp", sep=" ")))
dev.off()
}

