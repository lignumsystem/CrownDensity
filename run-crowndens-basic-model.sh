#NOTE: the dummy '-iter <YEARS>' and -metafile MetaFile.txt argument. HDF5 files are generated with GrowthLoop from LignumForest.
#The GrowthLoop command line checks for -iter argument.
#Otherwise `crowndens` has compatible command line with `lignum-forest`.
#NOTE: Glob expressions can be uses for MetaFiles: crowndens 60 '{Metafile,Metafile1}.txt'
#      matches MetaFile.txt and MetaFile1.txt and nothing else.
#      Protect the Glob expression with hyphens ('). See example in the command line below.
#NOTE: -modeChange accepts comma separated list of years, e.g. -modeChage 10,20,30
#NOTE: Consistency check is made: Number of MetaFiles must be number of ModeChange years + 1
#      (The first Metafile is consumed in the first initialization).
#REMEMBER: Set -hdf5 <filename> !!!!!!!
./crowndens 80 '{MetaFile,MetaFile1}.txt' -metafile MetaFile.txt -iter 80  -voxelspace VoxelSpace.txt -voxelCalculation 5 -numParts 2 -hw 5 -kBorderConifer 0.14 -modeChange 10 -architectureChange 10 -aChangeStart 1 -Lmaxturn 80 -writeInterval 5 -increaseXi 15 -dumpSelf -n_buds_ini_min 3 -n_buds_ini_max 3 -verbose -butt_swell_coeff 0.0015 -butt_swell_start 0 -hdf5 koe21.h5
