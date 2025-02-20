##N-body Simulation in C++

This program simulates the motion of particles under gravitational forces
, based on the laws of physics. The simulation computes the interaction
between particles, updates their velocities and positions, and outputs the 
state of the system at regular intervals.

Things to Have:
C++ Compiler 
Standard C++ library
Text editor or IDE

Step 1: Open Terminal (Command Prompt)
On Windows, use the command prompt or a terminal emulator like Git Bash.

Step 2: Navigate to the Project Directory
Use the cd command to navigate to the directory containing your code.

cd /path/to/nbody_simulation

Step 3: Compile the Code
To compile the code, run the following command (ensure that your compiler supports C++11 or later):
On Windows run this command:
g++ N-body-simulation.cpp -o N-body-simulation

This will produce an executable file named N-body-simulation

Step 4: Verify the Compilation
After the compilation process, check if the executable was created in the current directory by listing the files:
ls

You should see N-body-simulation

Step 5:Running the Simulation
To run the simulation, use the following command:

./N-body-simulation <num_particles> <dt> <num_steps> <dump_interval>

Where:
<num_particles>: Number of particles in the simulation (integer).
<dt>: Time step for the simulation (double).
<num_steps>: Number of time steps to run the simulation (integer).
<dump_interval>: The interval (in steps) to output the state of the simulation (double).
  
Example Command:
./nbody_simulation 100 0.01 1000 10

This runs the simulation with 100 particles, a time step of 0.01, for 1000 steps, and outputs the state of the simulation every 10 steps.
