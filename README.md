# hamFoam
Heat and Moisture transport solver for building physics. OpenFOAM based.

Tested on OpenFOAM-dev. Should work atleast on OpenFOAM-7 also.


Initial testing has been done. Works. Finer details about material properties etc. will be adjusted. Once the adjustements have been done, a comparison to COMSOL will be added.

Measurement data by Tampere University, Building Physics team led by Juha Vinha. We have more data than is published here.


# To compile
cd solvers
./Allwmake

# To test transient hamFoam
(for all the scripts to work, you need python3 and some packages)
cd validationAndVerification/hamFoam/
bash runAll.sh

Antti Mikkonen 2019, a.mikkonen@iki.fi
