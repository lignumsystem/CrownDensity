
#########        CompareTwoRuns function prints out graphs from two runs in HDF5 files
#########        Call: CompareTwoRuns(infile1, infile2)
#########        infile1,1: names of input HDF5 file (if not in the current directory include path). The file
#########                must be in the HDF5 File Format (see below)


#The Lignum simulation of one tree (in CrownDensity) is stored in a file in the HDF5 File Format
#(https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

#In order to run CompareTwoRuns in R a suitable package is needed. CompareTwoRuns has been tested using
#Biomanager package (https://www.rdocumentation.org/packages/BiocManager/versions/1.30.17):
#install.packages("BiocManager")
#BiocManager::install("rhdf5")

#To run CompareTwoRuns load the library first
#library("rhdf5")

#Then source "CompareTwoRuns.R"
#source("CompareTwoRuns.R")
#Finally run the CompareTwoRuns function, e.g.:
#CompareTwoRuns("CrownDensity147.h5","CrownDensity148.h5)
#The figures will appear in "koe.pdf"


CompareTwoRuns <- function(infile1, infile2) {

d1 <- H5Fopen(infile1)
d2 <- H5Fopen(infile2)


pdf("koe.pdf")

print(paste(infile1," = green",sep=""))
print(paste(infile2," = blue",sep=""))

td1 <- d1$ForestTreeData           #Tree Data
td2 <- d2$ForestTreeData           #Tree Data

ymax1 <- length(td1[1,1,]) - 1     #Maximum year
y1 <- 0:ymax1                      #Simulation years
ymax2 <- length(td2[1,1,]) - 1     #Maximum year
y2 <- 0:ymax2                      #Simulation years


#For stand
st1 = h5read(d1,"AllFunctions/stemsha.fun")
st2 = h5read(d2,"AllFunctions/stemsha.fun")

#d$AllFunctions$stemsha.fun.       #stemsha.fun function that speciefies the stand density

dens1 <- approx(st1[1,],st1[2,],y1)$y
dens2 <- approx(st2[1,],st2[2,],y2)$y

#Height
plot(y1,td1[7,1,], type="l", ylim=c(0,30), lwd=2, xlab="time (y)", ylab="Tree height (m)", main=paste("Tree height\n",infile1," = green  ",infile2," = blue",sep=""),col="darkgreen")
points(y2,td2[7,1,], type="l", lwd=2, lty=1,col="blue")
points(y1,td1[11,1,], type="l", lwd=2,col="darkgreen",lty=2)
points(y2,td2[11,1,], type="l", lwd=2, lty=2,col="blue")

va27 <- read.table("ResultAnalysis/Va27.txt",header=FALSE)
colnames(va27) <- c("a", "Hd",   "HgM",    "DgM",   "V",  "Hc", "G")
vv <- read.table("ResultAnalysis/VVV-40.txt",header=FALSE)
colnames(vv) <- c("age",    "DBH",              "H",            "Hcb",          "Wf", "V")
ksto <- read.table("ResultAnalysis/ksto-mt.dat", header=TRUE)
points(va27$a,va27$HgM,type="l",lwd=2,col="brown")
points(vv$age,vv$H,type="l",lwd=2,col="brown")
points(ksto$year,ksto$Hav,type="l",lwd=2,col="brown")
legend(1,30,"brown = Varmola, Vuokila&V:aho, Koivisto",box.lty=0,text.col="brown")



#Base diameter
plot(y1,100*td1[8,1,], type="l", lwd=2, ylim=c(0,30),xlab="time (y)", ylab="Base diameter (cm)",
main="Diameter at base",col="darkgreen")
points(y2,100*td2[8,1,], type="l", lwd=2,col="blue")

#Self thinning plot
plot(log(td1[8,1,]),log(dens1), xlim=c(log(0.001),log(0.7)),ylim=c(log(100),log(20000)),type="l", lty=1, lwd=2, xlab="log(base diameter (m))", ylab="log(No. trees / ha)", main="Self-thinning curve\nbrown = Koivisto", col = "darkgreen")
points(log(td2[8,1,]),log(dens2), type="l", lty=1, lwd=2, col="blue")
p1 <- c(max(log(td1[8,1,ymax1]))+1,min(log(dens1[ymax1]))-0.5)
p22 <- log(dens1[1])+0.5
p21 <- (p22-p1[2])/(-3/2)+p1[1]
points(c(p1[1],p21),c(p1[2],p22),type="l",lwd=2,col="red")
legend(p1[1]-0.5,p1[2]+1,"-3/2",box.lty=0,text.col="red")
points(log((0.02+1.25*ksto$Dbhav)/100),log(ksto$N),type="l",lwd=2,col="brown")


#Height vs diameter
plot(100*td1[8,1,],td1[7,1,], ylim=c(0,30), xlim=c(0,30), type="l",main="Height vs diameter at base" ,xlab="Tree diameter (m)", ylab="Tree height (cm)",col="darkgreen",lwd=2)
points(100*td2[8,1,],td2[7,1,], type="l", col="blue", lwd=2)
abline(0,1,lwd=2)

#Basal area of one tree
plot(y1,pi*(100*td1[8,1,])^2/4, ylim=c(0,500),type="l", lwd=2,xlab="time (y)", ylab= "cm2", main="Basal area at base of one tree",col="darkgreen")
points(y2,pi*(100*td2[8,1,])^2/4, type="l", lwd=2,col="blue")

#Basal area at crown base
plot(y1,pi*(100*td1[10,1,])^2/4, ylim=c(0,500),type="l", lwd=2,xlab="time (y)", ylab= "cm2", main="Basal area at crown base of one tree",col="darkgreen")
points(y2,pi*(100*td2[10,1,])^2/4, type="l", lwd=2,col="blue")

#Stem volume
plot(y1,td1[16,1,], ylim=c(0,1),type="l", lwd=2,xlab="time (y)", ylab= "m3", main="Stem volume of one tree",col="darkgreen")
points(y2,td2[16,1,], type="l", lwd=2,col="blue")


#Leaf area
plot(y1,td1[15,1,], ylim=c(0,300),type="l", lty=1, lwd=2,xlab="time (y)", ylab="m2",main="All-sided needle area of one tree",col="darkgreen") 
points(y2,td2[15,1,], type="l", lwd=2,col="blue")

#Photosynthetic production and respiration of one tree
plot(y1,td1[17,1,],type="l", lwd=2,xlab="time (y)", ylab= "kg C", main="Photos. (solid), Resp. (dashed)",ylim=c(0,10),col="darkgreen") #   P/Af
points(y2,td2[17,1,], type="l", lwd=2,col="blue")

points(y1,td1[18,1,], type="l", lty=2, lwd=2,col="darkgreen")
points(y2,td2[18,1,], type="l", col="blue", lwd=2,lty=2)

#Crown ratio
plot(y1,1-td1[11,1,]/td1[7,1,], ylim=c(0,1), type="l", main="Crown ratio", ylab="Crown ratio",xlab="time (y)",xlim=c(0,ymax1),col="darkgreen",lwd=2)
points(y2,1-td2[11,1,]/td2[7,1,], type="l", lwd=2,col="blue")


#Foliage mass vs cross-sectional area at crown base
plot((pi/4)*100^2*td1[10,1,]^2,2*td1[37,1,], ylim=c(0,10), xlim=c(0,150), type="l",main="Fol.mass vs stem cross sec. area",xlab="Stem cross section area at crown b.  (cm2)", ylab="Foliage mass (kg DM)",col="darkgreen",lwd=2)
points((pi/4)*100^2*td2[10,1,]^2,2*td2[37,1,],type="l", col="blue", lwd=2)
abline(0,0.055,lwd=2)


plot(y1,td1[35,1,]/(td1[37,1,]*td1[32,1,]), type="l",main="Radiation capture efficiency, related to STAR of tree", xlab="time (y)", ylab="Qabs/(Af*QinTop)",ylim=c(0,2),col="darkgreen",lwd=2)
points(y2,td2[35,1,]/(td2[37,1,]*td2[32,1,]), type="l", lwd=2,col="blue")


#========================================================================


#Basal area
plot(y1,pi*(td1[8,1,])^2*dens1/4, ylim=c(0,80),type="l", lwd=2,xlab="time (y)", ylab= "m2/ha", main="Stand basal area",col="darkgreen")
points(y2,pi*(td2[8,1,])^2*dens2/4,type="l", lwd=2, col="blue")

#Basal area at crown base
plot(y1,pi*(td1[10,1,])^2*dens1/4, ylim=c(0,60),type="l", lwd=2,xlab="time (y)", ylab= "m2/ha", main="Stand basal area at crown base",col="darkgreen")
points(y2,pi*(td2[10,1,])^2*dens2/4,type="l", lwd=2, col="blue")

#Stem volume
plot(y1,td1[16,1,]*dens1, ylim=c(0,1000),type="l", lwd=2,xlab="time (y)", ylab= "m3/ha", main="Stand stem volume",col="darkgreen")
points(y2,td2[16,1,]*dens2,type="l", lwd=2, col="blue")

#LAI and specific leaf area
#par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(y1,td1[15,1,]*dens1/1e4, ylim=c(0,30),type="l", lty=1, lwd=2,xlab="time (y)", ylab="All-sided LAI (m2/m2)",main="LAI (= Af * density)",col="darkgreen") #LAI
points(y2,td2[15,1,]*dens2/1e4,type="l",lty=1,lwd=2,col="blue")

#Photosynthetic production and respiration
plot(y1,td1[17,1,]*dens1/1000,type="l", lwd=2,xlab="time (y)", ylab= "t C / ha", main="Stand GPP (solid), Respiration (dashed)",ylim=c(0,15),col="darkgreen")
points(y2,td2[17,1,]*dens2/1000,type="l", lwd=2,col="blue")

points(y1,td1[18,1,]*dens1/1000, type="l", lty=2, lwd=2,col="darkgreen") 
points(y2,td2[18,1,]*dens2/1000, type="l", col="blue", lwd=2, lty=2)   #min




dev.off()

}
