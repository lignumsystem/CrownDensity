#NOTE: the dummy '-iter <YEARS>' and -metafile MetaFile.txt argument. HDF5 files are generated with GrowthLoop from LignumForest.
#The GrowthLoop command line checks for -iter argument.
#Otherwise `crowndens` has compatible command line with `lignum-forest`.
#REMEMBER: Set -hdf5 <filename> !!!!!!!
./crowndens 60 MetaFile.txt -metafile MetaFile.txt -iter 60  -voxelspace VoxelSpace.txt -voxelCalculation 5 -numParts 2 -hw 10 -kBorderConifer 0.11 -modeChange 20 -architectureChange 5 -aChangeStart 20 -writeInterval 5 -increaseXi 15  -hdf5 fip1-2fgo1-x1ac5.h5
#-modeChange 20 -architectureChange 5 -aChangeStart 20
#Ajo147
#clargs <- c("dummy", "34", "MetaFile.txt", "-numParts", "2", "-hw",  "10", "-voxelCalculation", "5", "-budViewFunction")
