#Rscript to generate single gene graph for ICR96 data
#Usage: Rscript PlotICR96single.r <id> <gene>

args = commandArgs(trailingOnly=TRUE)
id = args[1]
gene = args[2]

dummy1 <- read.csv(paste("Fx-",id,".csv",sep=""))
gpos = c(" ", 0, 0, " ")
if ( gene == "SDHB") {
gpos = c("chr1", 17018880, 17054020, "SDHB")
} else if ( gene == "MUTYH") {
gpos = c("chr1", 45329305, 45340255, "MUTYH")
} else if ( gene == "FH") {
gpos = c("chr1", 241497827, 241519723, "FH")
} else if ( gene == "EPCAM" ) {
gpos = c("chr2", 47369505, 47386614, "EPCAM")
} else if ( gene == "MSH2" ) {
gpos = c("chr2", 47403191, 47482950, "MSH2")
} else if ( gene == "MSH6" ) {
gpos = c("chr2", 47783233, 47806861, "MSH6")
} else if ( gene == "BARD1" ) {
gpos = c("chr2", 214728664, 214809581, "BARD1")
} else if ( gene == "MLH1" ) {
gpos = c("chr3", 36993547, 37050654, "MLH1")
} else if ( gene == "BAP1") {
gpos = c("chr3", 52402287, 52409879, "BAP1")
} else if ( gene == "APC") {
gpos = c("chr5", 112737876, 112844127, "APC")
} else if ( gene == "NSD1" ) {
gpos = c("chr5", 177133913, 177295460, "NSD1")
} else if ( gene == "PMS2" ) {
gpos = c("chr7", 5973398, 6009020, "PMS2")
} else if ( gene == "EZH2" ) {
gpos = c("chr7", 148807645, 148847299, "EZH2")
} else if ( gene == "NBN" ) {
gpos = c("chr8", 89935581, 89984562, "NBN")
} else if ( gene == "CDKN2A" ) {
gpos = c("chr9", 21968228, 21994332, "CDKN2A")
} else if ( gene == "BMPR1A") {
gpos = c("chr10", 86756627, 86923720, "BMPR1A")
} else if ( gene == "PTEN") {
gpos = c("chr10", 87864469, 87965473, "PTEN")
} else if ( gene == "WT1" ) {
gpos = c("chr11", 32389057, 32435346, "WT1")
} else if ( gene == "ATM" ) {
gpos = c("chr11", 108222820, 108365509, "ATM")
} else if ( gene == "CDK4" ) {
gpos = c("chr12", 57748524, 57752533, "CDK4")
} else if ( gene == "BRCA2" ) {
gpos = c("chr13", 32315462, 32398771, "BRCA2")
} else if ( gene == "RB1" ) {
gpos = c("chr13", 48303912, 48480072, "RB1")
} else if ( gene == "PALB2") {
gpos = c("chr16", 23603458, 23641158, "PALB2")
} else if ( gene == "CDH1") {
gpos = c("chr16", 68737415, 68833500, "CDH1")
} else if ( gene == "TP53" ) {
gpos = c("chr17", 7669608, 7687550, "TP53")
} else if ( gene == "NF1" ) {
gpos = c("chr17", 31095309, 31374156, "NF1")
} else if ( gene == "RAD51D" ) {
gpos = c("chr17", 35100952, 35119614, "RAD51D")
} else if ( gene == "BRCA1" ) {
gpos = c("chr17", 43045677, 43125382, "BRCA1")
} else if ( gene == "RAD51C" ) {
gpos = c("chr17", 58692643, 58734223, "RAD51C")
} else if ( gene == "BRIP1" ) {
gpos = c("chr17", 61683295, 61863406, "BRIP1")
} else if ( gene == "SMAD4" ) {
gpos = c("chr18", 51030201, 51078468, "SMAD4")
} else if ( gene == "STK11" ) {
gpos = c("chr19", 1206913, 1228447, "STK11")
} else if ( gene == "CHEK2" ) {
gpos = c("chr22", 28687896, 28741856, "CHEK2")
}

dummy2 <- dummy1[dummy1$Chr==gpos[1] & dummy1$Position>=gpos[2] & dummy1$Position<=gpos[3], ]

pdf(file = paste("Rplot",id,gene,".pdf"))
par(mar=c(6,5,5,5)+0.1)
plot(dummy2$Position,dummy2$copynumber,main=paste(id,gene),ylab="Copy Number", xlab=paste(dummy2$Chr[1]," (hg38)"), ylim=c(0,10))
abline(h=0,col="black",lwd=2)
abline(h=1.0,col="red",lwd=2)
abline(h=2.0,col="green",lwd=2)
abline(h=3.0,col="blue",lwd=2)
par(new=TRUE)
## Plot the second plot and put Z-score axis scale on right
plot(dummy2$Position, -abs(dummy2$zscore), pch=5, cex=0.75, xlab="", ylab="", ylim=c(-10,0), 
    axes=FALSE, type="b", col="red")
mtext("-abs(Z-score)",side=4,col="red",line=3) 
axis(4, ylim=c(-10,0), col="red",col.axis="red",las=1)
abline(h=0,lty=4,col="black")
abline(h=-1.67,lty=3,col="black")
abline(h=-2.3,lty=2,col="black")
dev.off()


