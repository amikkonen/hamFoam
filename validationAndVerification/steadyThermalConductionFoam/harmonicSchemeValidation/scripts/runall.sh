bash scripts/cleanRun.sh
bash scripts/mesh.sh

rm -rf postProcessing*

# Harmonic
cp system/fvSchemes_harmonic system/fvSchemes
bash scripts/solve.sh
cp -rf postProcessing/ postProcessing_harmonic
bash scripts/cleanRun.sh

# Linear
cp system/fvSchemes_linear system/fvSchemes
bash scripts/solve.sh
cp -rf postProcessing/ postProcessing_linear
bash scripts/cleanRun.sh

python3 validate.py

