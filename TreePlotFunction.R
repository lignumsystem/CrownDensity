
#########        ForestPlot function prints out graphs from Lignum simulation of a forest stand
#########        Call: ForestPlot(infile, aplot, pick)
#########        infile: name of input HDF5 file (if not in the current directory include path). The file
#########                must be in the HDF5 File Format (see below)
#########        aplot   Area of the simulated plot in m2
#########        pick    Every pick tree will appear in the height and diameter growth graphs.
#########                (If there are e.g. 600 trees on plot graphs become a mess if all are plotted)

## This TreePlotFunction.R is modified from ForestPlotFunction.R (LignumForest project) for a single tree (CrownDensity project).
## In short all StandData is removed and Leaf Area Index (LAI) and Specific Leaf Area (SLA) for a single tree
## based on single tree foliage area, single tree average branch length and single tree foliage mass.
## Data is available straightforwardly from collected TreeForestData.
## ----------------------------------------------------------------------------------------------------------------------------
##  The Lignum forest simulation results are stored in a file in the HDF5 File Format
##  (https://en.wikipedia.org/wiki/Hierarchical_Data_Format).
##
##  ForestPlot assumes the input file contains the following datasets: 
##                /  ForestTreeData H5I_DATASET  FLOAT 51 x 599 x 21
##                /       StandData H5I_DATASET  FLOAT       19 x 21
##
##  In this example 599 is the numberof trees on plot and 21 is the number times (years) at which
##  results are stored
##
##  In order to run ForestPlot in R a suitable package is needed. ForestPlot has been tested using
##  Biomanager package (https://www.rdocumentation.org/packages/BiocManager/versions/1.30.17):
##  install.packages("BiocManager")
##  BiocManager::install("rhdf5")
##
##  For Lorenz curve (inequality) install the package "ineq"
##  install.packages("ineq")
##  
##  To run the TreePlot function load libraries first
##  library("rhdf5")
##  library("ineq")
##  Then source "TreePlotFunction.R"
##  source("TreePlotFunction.R")
##  Finally run the ForestPlot function, e.g.:
##  TreePlot("CrownDensity147.h5",1,1)
##  The figures will appear in "../CrownDensity147.h5.pdf"
##
##  NOTE! Three data files are are used: Va27.txt, VVV-40.txt, and ksto-mt.dat (in LignumForest/Resultanalysis).
##  Va27.txt = Varmola M (1987) Männyn viljelytaimikoiden kasvumalli. Lic. For & Agric.
##  Thesis, Department of Forest Mensuration, University of Helsinki, 89 p.
##  VVV-40.txt = Vuokila Y, Väliaho H (1980) Viljeltyjen havumetsiköiden kasvatusmallit.
##  Communicationes Instituti Forestalis Fenniae 99, 271.
##  ksto-mt.dat.  Finnish yiels tables compiled by Koivisto
##   -- You will need to adjust the path to these files below if you are not running this function
##   in LignumForest/Resultanalysis (getwd() == LignumForest/Resultanalysis)
##  ----------------------------------------------------------------------------
gini <- function(v){
    v <- na.omit(v)
    n <- length(v)
    sum(outer(v, v, FUN=function(x,y){abs(x-y)}),na.rm=TRUE) / (2 * n * sum(v,na.rm=TRUE))
}

TreePlot <- function(infile,aplot,pick) {

    ## Open infile and create outfile 
    d <- H5Fopen(infile)

    pdf_file <- paste(infile,".pdf",sep="")
    pdf(pdf_file)

    ## For the 'y' take the number of years simulated
    dset=d$ForestTreeData
    dims = dim(dset)
    simyears <- dims[3]
    y = seq(from=1,simyears)
    ## Collect longest, tallest median trees. Not quite meaningful in this single tree context
    ## but used in plotting longest and shortest trees
    h <- d$ForestTreeData[7,,]
    ymax <- max(h,na.rm=TRUE)
    mh <- max(d$ForestTreeData[7,,ymax],na.rm=TRUE)
    largest <- which(h>0.999*mh)[1]
    print(largest)
    mh <- min(d$ForestTreeData[7,,ymax],na.rm=TRUE)
    smallest <- which(h<1.001*mh)[1]
    
    mh <- median(d$ForestTreeData[7,,ymax],na.rm=TRUE)
    med <- which(h<1.002*mh&h>0.98*mh)[1]

    ##Plotting starts
    ##Height
    print("Height")
    plot(y,d$ForestTreeData[7,1,], type="l",xlab="time (y)", ylab="Tree height (m)",lwd=2,
         main=paste("Height growth"))
    print("Base diameter")
    ##Base diameter
    plot(y,100*d$ForestTreeData[8,1,], type="l",ylab="Diameter at base (cm)",lwd=2,
         main=paste("Base diameter growth"))

    print("Max branch extension")
    plot(y,d$ForestTreeData[52,1,], type="l",ylab="Max branch extension (m)",lwd=2,
         main=paste("Branch extension"))
    print("LAI with max branch")
    ##LAI and specific leaf area for a single tree
    par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
    ##LAI in Forest 
    ##plot(y,d$StandData[17,], ylim=c(0,20),type="l", lty=1, lwd=2,xlab="time (y)", ylab="All-sided LAI (m2/m2)",main="LAI = continuous,  specific LA = dashed")
    ##LAI for single tree 15="TreeAf", 50 = "MeanBranch_SumL/Nbranch", meant to be a mean branch that is horizontal 
    plot(y,d$ForestTreeData[15,1,]/(pi*(d$ForestTreeData[52,1,]**2.0)), type="l", lty=1, lwd=2, 
         xlab="time (y)", ylab="LAI (m2/m2) as TreeAf/pi*MaxBranchExtension^2",main="LAI with max branch extension = continuous,  specific LA = dashed")
    par(new = TRUE)
    ##SLA in Forest
    ##plot(y,d$StandData[17,]/d$StandData[18,],type="l", lty=2,lwd=2, axes=FALSE,xlab = "", ylab = "", ylim=c(0,32))
    ##SLA for single tree, 15 = "TreeAf", 23="Wf"
    print("SLA")
    plot(y,d$ForestTreeData[15,1,]/d$ForestTreeData[23,1,],type="l", lty=2,lwd=2, axes=FALSE,xlab = "", ylab = "", ylim=c(0,32))
    axis(side = 4, at = pretty(c(0,32)))
    mtext("SLA (m2/ kg C) as TreeAf/TreeWf", side = 4, line = 3)             # Add second axis label
    print("LAI with mean branch")
    plot(y,d$ForestTreeData[15,1,]/(pi*(d$ForestTreeData[50,1,]**2.0)), type="l", lty=1, lwd=2,
         xlab="time (y)", ylab="LAI (m2/m2) as TreeAf/pi*MeanBranch_SumL/Nbranch^2",main="LAI with mean branch = continuous,  specific LA = dashed")
    par(new = TRUE)
    print("SLA")
    ##SLA in Forest
    ##plot(y,d$StandData[17,]/d$StandData[18,],type="l", lty=2,lwd=2, axes=FALSE,xlab = "", ylab = "", ylim=c(0,32))
    ##SLA for single tree, 15 = "TreeAf", 23="Wf"
    plot(y,d$ForestTreeData[15,1,]/d$ForestTreeData[23,1,],type="l", lty=2,lwd=2, axes=FALSE,xlab = "", ylab = "", ylim=c(0,32))
    axis(side = 4, at = pretty(c(0,32)))
    mtext("SLA (m2/ kg C) as TreeAf/TreeWf", side = 4, line = 3)             # Add second axis label
    
    print("P and M")
    ##Photosynthetic production and respiration single tree
    ##P in forest
    ##plot(y,apply(d$ForestTreeData[17,,],2,sum,na.rm=TRUE),type="l", lwd=2,xlab="time (y)", ylab= "kgC", main="GPP and Respiration") #
    ##P for single tree
    plot(y,d$ForestTreeData[17,1,],type="l", lwd=2,xlab="time (y)", ylab= "kgC", main="GPP and Respiration")
    ##M for forest
    ##points(y,apply(d$ForestTreeData[18,,],2,sum,na.rm=TRUE), type="l", lty=2, lwd=2)   #min
    ##M for single tree
    points(y,d$ForestTreeData[18,1,],type="l", lty=2, lwd=2)
    
    print("TreeH and TreeHCrownBase")
    ##Tree heights and Crown heights
    plot(y,d$ForestTreeData[7,1,],type="l", main=paste("Tree and Tree crown heights"),
         ylab="Tree height (m)",xlab="time (y)")
    points(y,d$ForestTreeData[11,1,],type="l",lty=2)

    print("TreeDbase")
    ##Tree diameters at base
    plot(y,100*d$ForestTreeData[8,1,], ylim=c(0,30), type="l",
         main="Individual tree diameter at base",xlab="time (y)", ylab="Tree diameter (cm)")

    print("Crown ratio")
    ##Crown ratio
    plot(y,1-d$ForestTreeData[11,1,]/d$ForestTreeData[7,1,], ylim=c(0,1), type="l",
         main="Crown ratios and mean", ylab="Crown ratio",xlab="time (y)",xlim=c(0,ymax))

    print("Height vs diameter")
    ##Height vs diameter
    plot(100*d$ForestTreeData[8,1,],d$ForestTreeData[7,1,], ylim=c(0,30), type="l",main="Height vs diameter at base",
         xlab="Tree diameter (m)", ylab="Tree height (cm)")
    abline(0,1,lwd=2)

    print("Foliage mass")
    plot(y,d$ForestTreeData[23,1,],type="l",
         main="Foliage mass",
         xlab="time(y)", ylab="Foliage mass (kg DM)")
    print("Diameter Crown base")
    plot(y,d$ForestTreeData[10,1,],type="l",main="Diameter crown base",
         ylab="Diameter crown base (m)",xlab="time (y)")
    print("Cross sectional area crown base")
    plot(y,pi*(d$ForestTreeData[10,1,]/2.0)**2,type="l",main="Cross sectional area crown base",
         ylab="Cross section area (m^2)",xlab="time (y)") 
    print("Foliage mass vs Cross sectional area crown base")
    ##Foliage mass vs cross-sectional area at crown base
    plot(100.0**2*pi*(d$ForestTreeData[10,1,]/2.0)**2,2.0*d$ForestTreeData[23,1,],type="l",
         main="Foliage mass vs stem cross sectional area at crown base",
         xlab="Stem cross section area at crown base.  (cm2)", ylab="Foliage mass (kg DM)")
    abline(0,0.055,col="blue",lwd=2)

    print("P/Af")
    ##P / Af
    plot(y,d$ForestTreeData[17,1,]/d$ForestTreeData[15,1,],type="l", lwd=2,xlab="time (y)",
         ylab= "P/Af min, mean, max (kg C / MJ PAR)", main="Photosynthetic rate per unit needle area") #   P/Af
  
    print("Radiation capture efficency as Qabs/TreeAf")
    ##Radiation capture efficiency
    plot(y,d$ForestTreeData[35,1,]/(d$ForestTreeData[32,1,]*d$ForestTreeData[15,1,]), type="l",
         main="Radiation capture efficiency, related to STAR of tree", xlab="time (y)", ylab="Qabs/(Af*QinTop)",ylim=c(0,0.2))
    points(y,d$ForestTreeData[35,1,]/(d$ForestTreeData[32,1,]*d$ForestTreeData[15,1,]), type="l")
    
    print("Plotting done")
 
    dev.off()

}
