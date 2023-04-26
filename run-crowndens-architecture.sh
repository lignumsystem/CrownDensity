#This is used to run architecture effects. Note the missing arguments for EBH, adhoc, random etc options. 
#NOTE The GrowthLoop from LignumForest is only used to crete HDF5 files.
#GrowthLoop checks for required command line arguments needed in LignumForest (acquired by `crowndens` in other ways)
#The GrowthLoop command line checks for -iter, -metafile and -kBorderConifer arguments. 
#HDF5 files need -hdf5 and -writeInterval command line arguments
#`crowndens` checks the two first positional arguments as number of iterations (years) and MetaFile respectively.
#Otherwise `crowndens` has compatible command line with `lig-forest`. 
./crowndens 60 MetaFile.txt -metafile MetaFile.txt -iter 60 -voxelspace VoxelSpace.txt  -voxelCalculation 5 -numParts 2 -hw 10 -kBorderConifer 0.11  -verbose  -writeInterval 5   -hdf5 CrownDensityArchitecture170423_20_147.h5
#-Ebhreduction 0.95 -EBHFINAL 0.5 -EBHInput 2
#Ajo147
#clargs <- c("dummy", "34", "MetaFile.txt", "-numParts", "2", "-hw",  "10", "-voxelCalculation", "5", "-budViewFunction")
