
#########        TreePlot function prints out graphs from Lignum simulation of one tree
#########		 (analogous to ForestPlot in LignumForest)
#########        Call: TreePlot(infile)
#########        infile: name of input HDF5 file (if not in the current directory include path). The file
#########                must be in the HDF5 File Format (see below)


#The Lignum simulation of one tree (in CrownDensity) is stored in a file in the HDF5 File Format
#(https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

#ForestPlot assumes the input file contains the following dataset: 
#              /  ForestTreeData H5I_DATASET  FLOAT 51 x 1 x simulation years

#In order to run ForestPlot in R a suitable package is needed. ForestPlot has been tested using
#Biomanager package (https://www.rdocumentation.org/packages/BiocManager/versions/1.30.17):
#install.packages("BiocManager")
#BiocManager::install("rhdf5")

#To run TreePlot load the library first
#library("rhdf5")

#Then source "TreePlotFunction.R"
#source("TreePlotFunction.R")
#Finally run the TreePlot function, e.g.:
#TreePlot("CrownDensity147.h5")
#The figures will appear in "CrownDensity147.h5.pdf"

#NOTE! Three data files are are used: Va27.txt, VVV-40.txt, and ksto-mt.dat (in LignumForest/Resultanalysis).
#Va27.txt = Varmola M (1987) Männyn viljelytaimikoiden kasvumalli. Lic. For & Agric.
#Thesis, Department of Forest Mensuration, University of Helsinki, 89 p.
#VVV-40.txt = Vuokila Y, Väliaho H (1980) Viljeltyjen havumetsiköiden kasvatusmallit.
#Communicationes Instituti Forestalis Fenniae 99, 271.
#ksto-mt.dat.  Finnish yiels tables compiled by Koivisto
# -- You will need to adjust the path to these files below if you are not running this function
# in LignumForest/Resultanalysis (getwd() == LignumForest/Resultanalysis)
#----------------------------------------------------------------------------

TreePlot <- function(infile) {

d <- H5Fopen(infile)

# va27 <- read.table("Va27.txt",header=FALSE)
# colnames(va27) <- c("a", "Hd",   "HgM",    "DgM",   "V",  "Hc", "G")

# vv <- read.table("VVV-40.txt",header=FALSE)
# colnames(vv) <- c("age",    "DBH",		"H",		"Hcb",		"Wf", "V")

# ksto <- read.table("ksto-mt.dat", header=TRUE)

pdf_file <- paste(infile,".pdf",sep="")
pdf(pdf_file)

td <- d$ForestTreeData           #Tree Data
ymax <- length(td[1,1,]) - 1     #Maximum year
y <- 0:ymax                      #Simulation years

#Height
plot(y,td[7,1,], type="l", ylim=c(0,30), lwd=2, xlab="time (y)", ylab="Tree height (m)", main="Tree height")
points(y,td[11,1,], type="l", lwd=2, lty=2)

# points(va27$a,va27$HgM,type="l",lwd=3,col="darkgreen")
# points(vv$age,vv$H,type="l",lwd=3,col="darkgreen")
# points(ksto$year,ksto$Hav,type="l",lwd=3,col="darkgreen")


#Base diameter
plot(y,100*td[8,1,], type="l", lwd=2, ylim=c(0,30),xlab="time (y)", ylab="Base diameter (cm)",
main="Diameter at base (black), breast height (brown), crown base (green)")
points(y,100*td[9,1,], type="l", lwd=2, col="brown")         #Diameter at breast height
points(y,100*td[10,1,], type="l", lwd=2, col="darkgreen")    #Diameter at crown base

# dkanto = 2 + 1,25d (Laasasenaho 1975, Folia Forestalia 233)
# points(va27$a,0.02+1.25*va27$DgM,type="l",lwd=3,col="darkgreen")
# points(vv$age,0.02+1.25*vv$DBH,type="l",lwd=3,col="darkgreen")
# points(ksto$year,0.02+1.25*ksto$Dbhav,type="l",lwd=3,col="darkgreen")


#Basal area of one tree
plot(y,pi*(100*td[8,1,])^2/4, ylim=c(0,500),type="l", lwd=2,xlab="time (y)", ylab= "cm2", main="Basal area at base of one tree")


#Basal area at crown base
plot(y,pi*(100*td[10,1,])^2/4, ylim=c(0,500),type="l", lwd=2,xlab="time (y)", ylab= "cm2", main="Basal area at crown base of one tree")

#Stem volume
plot(y,td[16,1,], ylim=c(0,1),type="l", lwd=2,xlab="time (y)", ylab= "m3", main="Stem volume of one tree")


#Leaf area
plot(y,td[15,1,], ylim=c(0,300),type="l", lty=1, lwd=2,xlab="time (y)", ylab="m2",main="All-sided needle area of one tree") 

#Photosynthetic production and respiration of one tree
plot(y,td[17,1,],type="l", lwd=2,xlab="time (y)", ylab= "kg C", main="Photos. (solid), Resp. (dashed), P - R (green) of one tree",
ylim=c(0,10)) #   P/Af
points(y,td[18,1,], type="l", lty=2, lwd=2)
points(y,td[17,1,]-td[18,1,], type="l", col="darkgreen", lwd=2)

#"GPP (solid), Respiration (dashed), GPP - Resp (green)"

#Crown ratio
plot(y,1-td[11,1,]/td[7,1,], ylim=c(0,1), type="l", main="Crown ratio", ylab="Crown ratio",xlab="time (y)",xlim=c(0,ymax))

#Height vs diameter
plot(100*td[8,1,],td[7,1,], ylim=c(0,30), xlim=c(0,30), type="l",main="Height vs diameter at base" ,xlab="Tree diameter (m)", ylab="Tree height (cm)")
abline(0,1,lwd=2)

#Foliage mass vs cross-sectional area at crown base
plot((pi/4)*100^2*td[10,1,]^2,2*td[37,1,], ylim=c(0,10), xlim=c(0,150), type="l",main="Fol.mass vs stem cross sec. area",xlab="Stem cross section area at crown b.  (cm2)", ylab="Foliage mass (kg DM)")
abline(0,0.055,col="blue",lwd=2)


#P / Af
plot(y,td[17,1,]/td[37,1,], ylim=c(0,3), type="l", lwd=2,xlab="time (y)", ylab= "P/Af min, mean, max (kg C / m2)", main="Photosynthetic rate per unit needle area") #   P/Af


plot(y,td[35,1,]/(td[37,1,]*td[32,1,]), type="l",main="Radiation capture efficiency, related to STAR of tree", xlab="time (y)", ylab="Qabs/(Af*QinTop)",ylim=c(0,2))

#Lambda
plot(y,td[51,1,], type="l", lty=2, lwd=2,xlab="time (y)", ylab= expression(lambda),ylim=c(0,7), main = expression(paste(lambda," in Eq: New growth(",lambda,") = P - M")))

#========================================================================
#For stand
st = h5read(d,"AllFunctions/stemsha.fun")
#d$AllFunctions$stemsha.fun.       #stemsha.fun function that speciefies the stand density

dens <- approx(st[1,],st[2,],y)$y
plot(y,dens, type="l", ylim=c(0,20000),lty=1, lwd=2, xlab="time (y)", ylab="No. trees / ha", main="Stand density")

#Basal area
plot(y,pi*(td[8,1,])^2*dens/4, ylim=c(0,80),type="l", lwd=2,xlab="time (y)", ylab= "m2/ha", main="Basal area")
# points(va27$a,va27$G,type="l",lwd=3,col="darkgreen")
# points(ksto$year,ksto$G,type="l",lwd=3,col="darkgreen")

#Basal area at crown base
plot(y,pi*(td[10,1,])^2*dens/4, ylim=c(0,60),type="l", lwd=2,xlab="time (y)", ylab= "m2/ha", main="Basal area at crown base")

#Stem volume
plot(y,td[16,1,]*dens, ylim=c(0,1000),type="l", lwd=2,xlab="time (y)", ylab= "m3/ha", main="Stem volume")
# points(va27$a,va27$V,type="l",lwd=3,col="darkgreen")
# points(vv$age,vv$V,type="l",lwd=3,col="darkgreen")
# points(ksto$year,ksto$V,type="l",lwd=3,col="darkgreen")

#LAI and specific leaf area
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(y,td[15,1,]*dens/1e4, ylim=c(0,80),type="l", lty=1, lwd=2,xlab="time (y)", ylab="All-sided LAI (m2/m2)",main="LAI = continuous,  specific LA = dashed") #LAI
par(new = TRUE)
plot(y,td[15,1,]/td[37,1,],type="l", lty=2,lwd=2, axes=FALSE,xlab = "", ylab = "", ylim=c(0,32))
axis(side = 4, at = pretty(c(0,32)))
mtext("Specific leaf area (m2/ kg C)", side = 4, line = 3)             # Add second axis label

#Photosynthetic production and respiration
plot(y,td[17,1,]*dens/1000,type="l", lwd=2,xlab="time (y)", ylab= "t C / ha", main="GPP (solid), Respiration (dashed), GPP - Resp (green)",ylim=c(0,15)) #   P/Af
points(y,td[18,1,]*dens/1000, type="l", lty=2, lwd=2)   #min
points(y,td[17,1,]*dens/1000-td[18,1,]*dens/1000, type="l", col="darkgreen", lwd=2)   #min



dev.off()

}
