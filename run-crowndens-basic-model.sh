#NOTE: the dummy '-iter <YEARS>' and -metafile MetaFile.txt argument. HDF5 files are generated with GrowthLoop from LignumForest.
#The GrowthLoop command line checks for -iter argument.
#Otherwise `crowndens` has compatible command line with `lignum-forest`.
#REMEMBER: Set -hdf5 <filename> !!!!!!!
./crowndens 60 MetaFile.txt -metafile MetaFile.txt -iter 60 -voxelspace VoxelSpace.txt  -voxelCalculation 5 -modeChange 20 -architectureChange 10 -fipmode fip1.fun -fgomode fgo1.fun -numParts 2 -hw 10 -kBorderConifer 0.11  -writeInterval 5 -increaseXi 15  -hdf5 CrownDensity-147-Basic-Aug23-74.h5
#-Ebhreduction 0.95 -EBHFINAL 0.5 -EBHInput 2
#Ajo147
#clargs <- c("dummy", "34", "MetaFile.txt", "-numParts", "2", "-hw",  "10", "-voxelCalculation", "5", "-budViewFunction")
