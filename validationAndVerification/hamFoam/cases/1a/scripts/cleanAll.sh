rm -rf postProcessing/ processor* [1-9]* 0 log/solver; cp -rf 0.orig/ 0
foamCleanPolyMesh
rm -rf 0 log/*
rm *.foam
rm *.OpenFOAM
rm *.pdf

