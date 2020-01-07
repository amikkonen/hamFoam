# hamFoam
Heat and Moisture transport solver for building physics. OpenFOAM based.

Tested on OpenFOAM-dev. Should work atleast on OpenFOAM-7 also.

Initial testing has been done. Works. Finer details about material properties etc. will be adjusted. Once the adjustements have been done, a comparison to COMSOL will be added.

Measurement data by Tampere University, Building Physics team led by Prof. Juha Vinha. More data is available than is published here.

# Previous publications on the measurements:

Vinha, J. 2007. Hygrothermal Performance of Timber - Framed External Walls in Finnish Climatic Conditions: A Method for Determining the Sufficient Water Vapour Resistance of the Interior Lining of a Wall Assembly. Publication 658. Tampere, Tampere University of Technology. 348 p

Kalamees, T. & Vinha, J. 2003. Hygrothermal calculations and laboratory tests on timber - framed wall structures. Building and Environment, e-published 8 February 2003, Vol. 38 (5), pp. 689-697.

Laukkarinen, A. & Vinha, J. 2011. Comparison of calculated and measured values of wall assembly tests using Delphin 5. Proceedings of the 9th Nordic Symposium on Building Physics, NSB 2011, Tampere, Finland, May 29 - June 2, Vol. 1, pp. 155-162.

# To compile
cd solvers
./Allwmake

# To test transient hamFoam
(for all the scripts to work, you need python3 and some packages)
cd validationAndVerification/hamFoam/
bash runAll.sh

Antti Mikkonen 2019, a.mikkonen@iki.fi
