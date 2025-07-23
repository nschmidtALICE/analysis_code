# Build the executable
make clean && make

# Run the compiled program
./bin/D0RecoEfficiency

# With environment variables
#FORCE_SERIAL=1 ./bin/D0RecoEfficiency
OMP_NUM_THREADS=4 ./bin/D0RecoEfficiency