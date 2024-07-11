##############################################################################################################
######### This script plots distribution of GC % at the center of provided sequences in the window of provided
######### length and compares it to the distribution across the genome.
#########   
######### input - seq_file_name      path to the FASTA file of DNA sequences to be avaluated
#########         genome_file_name   path to file with genome sequence in the format 1 in substituted for C/G and
#########                            0 substituted for A/T
#########         window            number of bp of the window to calculate GC %
##############################################################################################################
                                         
library(zoo)
library(sm)
library(optparse)
library(stringr) 

option_list = list(
make_option(c("-g", "--genome_file"), type="character", default=NULL, help="file with genome sequence", metavar="character"),
                   make_option(c("-f", "--file"), type="character", default=NULL, help="file with sample sequences", metavar="character"),
                   make_option(c("-w", "--window"), type="integer", default=200, help="size of the window to compare GC%", metavar="integer")
 );

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$file) || is.null(opt$genome_file)){
  print_help(opt_parser)
  stop("Genome and sample sequences files must be supplied (input file).n", call.=FALSE)
}

genome_file_name<-opt$genome_file
seq_file_name<-opt$file
out_file_name<-paste(seq_file_name,"_bias.jpg",sep="")
window<-opt$window

{
my_fasta<-readLines(seq_file_name)

if(length(my_fasta) < 2){stop("Bad fasta.")}
sss<-my_fasta[seq(2,length(my_fasta),2)]

t<-readLines(genome_file_name)


############## Verify user input
if(window != round(window)){stop("Window should be integer.")}
if(length(which(str_count(sss,"A|a|T|t|C|c|G|g") != nchar(sss))) != 0){stop("Sequences contain characters other than A, T, C and G.")}

rr<-strsplit(t[1],"")
rr_n<-as.numeric(unlist(rr))

##Calculate GC% in a rolling window of provided size
TS<-zoo(rr_n)
pp<-rollapply(TS, width = window, by = window, FUN = mean, na.rm = TRUE, align = "left")
center_of_region<-round(nchar(sss[1])/2)

sss_bin <- gsub("(A|T)", "0", sss, ignore.case = TRUE)
sss_bin <- gsub("(C|G)", "1", sss_bin, ignore.case = TRUE)

df<-as.data.frame(strsplit(sss_bin,split=""))
y<-colMeans(sapply(df[(center_of_region-round(window/2)):(center_of_region+round(window/2)),],as.numeric))


label_1<-rep(1,length(pp))
label_2<-rep(2,length(y))
label<-c(label_1,label_2)
dff<-c(as.numeric(pp),as.numeric(y))

###### Calculate p- value
ttest<-t.test(pp,y)
pval<-ttest[["p.value"]]
if(pval == 0){pval<-"p-value < 2.2e-16"}

jpeg(out_file_name,width=300, height=300)
  sm.density.compare(dff,label,xlab="GC %")
  title(main="GC % of sample versus\n genomic")
  legend("topright", legend = c("genome","sample"),
       lty = c(1,2), col = c("red","green"))
  text(0.1,5.5,pval)
dev.off()
}


# Hello pull request!
