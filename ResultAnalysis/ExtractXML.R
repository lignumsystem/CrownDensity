#########        ExtractXML function reads Lignum a tree stored in a HDF5 file produced by CrownDensity
#########        The function stores the XML file of the tree at given times.
#########        The names of the XML files are xmfile_Tree0_year.xml
#########        Call: ForestPlot(xmlfile, years)
#########        xmlfile:  HDF5 file containing the XML data
#########        years     vector of times for storing the XML files.

#The datafile must be stored iin the HDF5 File Format
#(https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

#The xmlfile must contain the XML of the tree as follows:
#           group     name       otype dclass dim
#           /TreeXML/5  Tree_0 H5I_DATASET STRING   1
#This example shows the xml data at time (year) 5

#In order to run ExtractXML in R a suitable package is needed. ExtractXML has been tested using
#Biomanager package (https://www.rdocumentation.org/packages/BiocManager/versions/1.30.17):
#install.packages("BiocManager")
#BiocManager::install("rhdf5")
#Library rhdf5
#library("rhdf5")




#----------------------------------------------------------------------------

ExtractXML <- function(xmlfile, years) {


#Write the XML strings to files
x <- H5Fopen(xmlfile)

for(i in 1:length(years)) {
	t <- paste("TreeXML/",as.character(years[i]),"/Tree_",as.character(0),sep="")
	if(H5Lexists(x,t)) {
		tree <- h5read(x,t)
		oname <- paste(xmlfile,"_Tree_",as.character(years[i]),".xml",sep="")
		cat(tree,file=oname)
		} else {
			print(paste("Tree ",t," not in the XML file",sep=""))
		}
	}
}

