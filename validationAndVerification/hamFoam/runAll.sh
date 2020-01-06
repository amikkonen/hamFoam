# Generate cases
echo "Generate cases"
python3 generate_1d_cases.py

# Run all cases
for i in 1a 1b 2a 2b 3a 3b 4a 4b
do
    echo "Running: cases/$i"
    cd cases/$i
    bash scripts/runall.sh
    #tsp -N 1 -L $(pwd) bash scripts/runall.sh 
    cd ../..
done

# Post process. Results in cases.
echo "Post processing"
python3 1DPost.py
echo "Ready. Results in: cases/all.pdf"
