#NOTE the dummy '-iter 4' and -metafile MetaFile.txt argument. HDF5 files are generated with GrowthLoop from LignumForest.
#The GrowthLoop command line checks for -iter argument.
#Otherwise `crowndens` has compatible command line with `lig-forest`. 
./crowndens 60 MetaFile.txt -metafile MetaFile.txt -iter 60 -voxelspace VoxelSpace.txt  -voxelCalculation 5 -numParts 2 -hw 10 -kBorderConifer 0.11 -EBH -verbose  -writeInterval 5 -increaseXi 15 -adHoc -RUE 1.0 -EBHREDUCTION 0.99 -EBHFINAL 0.5 -EBHInput 2 -hdf5 CrownDensity147.h5 
#-Ebhreduction 0.95 -EBHFINAL 0.5 -EBHInput 2

#Ajo147
#clargs <- c("dummy", "34", "MetaFile.txt", "-numParts", "2", "-hw",  "10", "-voxelCalculation", "5", "-budViewFunction")
