#!/bin/bash
#SBATCH --job-name=N-body-simulation
#SBATCH --partition=Centaurus
#SBATCH --time=01:00:00
#SBATCH --mem=32GB
#SBATCH --output=simulation_output.log

echo "Starting simulations..."

# Run solar  simulation and save to sim_solar_simulation_time.txt
echo "Running simulation: solar simulation, dt=200, steps=5000000"
(time ./N-body-simulation 12 200 5000000 1000 > solar_simulation.tsv) 2>&1 | tee sim_solar_simulation_time.txt

# Run 100-particle simulation and save to sim_100_simulation_time.txt
echo "Running simulation: 100 particles, dt=1, steps=10000"
(time ./N-body-simulation 100 1 10000 100 > solar_100.tsv) 2>&1 | tee sim_100_simulation_time.txt

# Run 1000-particle simulation and save to sim_1000_simulation_time.txt
echo "Running simulation: 1000 particles, dt=1, steps=10000"
(time ./N-body-simulation 1000 1 10000 100 > solar_1000.tsv) 2>&1 | tee sim_1000_simulation_time.txt

echo "Simulations completed."


