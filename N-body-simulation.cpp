#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

/**
TODO. Implement the N-body simulation in C++ by completing the following tasks:
 1. Define the Simulation Structure: Create a structure or class to represent the state of the simula
tion. This should include:
 1
• Particle properties: masses, positions, velocities, and forces (e.g., arrays or vectors to store these
 values for all particles).
 • Initialization functions:– Random initialization of particle properties (masses, positions, velocities).– Predefined configurations such as a simple 2 or 3 particle setup (e.g., Sun, Earth, Moon).– Load from file. Check recommended format.
 2. Force Calculation: Write a function to compute the gravitational force between every pair of particles
 and sum the force vectors for each particle. Ensure that forces are reset to zero at the start of each
 time step.
 3. Integration of Motion: Implement functions to:
 • Update particle velocities based on the computed forces and their masses.
 • Update particle positions based on their velocities and the time step.
 4. Output State: Outputs the state of the simulation (e.g., positions, velocities, forces) to a log file or
 standard output. If you use the recommended format for input and output, you can use the python
 script to visualize the output.
 5. Simulation Loop: Write a main simulation loop that:
 • Iterates over a fixed number of time steps.
 • Updates forces, velocities, and positions at each step.
 • At regular intervals, output the state to a log file.
 */

//Initialization Functions
//Random Initialization
//Predetermined Configurations
//Load from File
//Force Calculation
//Integration of Motion
//Output State
//Simulation Loop

//Gravitational Force
const double G = 6.67 * pow(10, -11);
//Softening Factor
const double epsilon = .0001;


//Particle Properties
struct Particle {
    double mass;
    double x, y, z; // position
    double vx, vy, vz; // velocity
    double fx, fy, fz; // force
};

/**
 * Resets the forces on each particle to zero.
 * This function is used to initialize or reset the force accumulators for each particle
 * before calculating the new forces in each simulation step.
 * 
 * @param particles A vector of Particle objects representing all particles in the simulation.
 */
void resetForces(std::vector<Particle>& particles) {
    for (auto& particleForce : particles) {
        // Set the force components to zero
        particleForce.fx = 0.0;
        particleForce.fy = 0.0;
        particleForce.fz = 0.0;
    }
}


/**
 * Calculates the forces between all particles in the simulation.
 * This function is used to calculate the gravitational forces between all particles in the simulation
 * before updating the velocities and positions of the particles.
 * 
 * @param particles A vector of Particle objects representing all particles in the simulation.
 */
void CalculateForce(std::vector<Particle>& particle) {
    // Call the resetForces function to reset the forces
    resetForces(particle);
    
    for (int i = 0; i < particle.size(); i++) {
        for (int j = i + 1; j < particle.size(); j++) {
            // Calculate the distance between the two particles
            // distance = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
            double dx = particle[i].x - particle[j].x;
            double dy = particle[i].y - particle[j].y;
            double dz = particle[i].z - particle[j].z;
            // r squared = r2 = r * r
            // distance = r = sqrt(r2)
            // Softening Factor used to prevent divisions by zero
            double r2 = ((dx * dx) + (dy * dy) + (dz * dz)) + epsilon;
            double r = std::sqrt(r2);
            // Calculate the force
            // m1 = particle[i].mass
            // m2 = particle[j].mass
            // F = G * m1 * m2 / r^2
            double force = (G * particle[i].mass * particle[j].mass)/ r2;

            // Calculate the direction of the force f =  d / r
            double fx = dx / r;
            double fy = dy / r;
            double fz = dz / r;

            // Add the force to particle 1
            particle[i].fx = particle[i].fx + force * fx;
            particle[i].fy = particle[i].fy + force * fy;
            particle[i].fz = particle[i].fz + force * fz;
            // Add the force to particle 2
            particle[j].fx = particle[j].fx - force * fx;
            particle[j].fy = particle[j].fy - force * fy;
            particle[j].fz = particle[j].fz - force * fz;
        }
    }
}

/**
 * Updates the velocities of all particles in the simulation.
 * This function is used to update the velocities of all particles in the simulation
 * using the calculated forces from the previous step.
 * 
 * @param particles A vector of Particle objects representing all particles in the simulation.
 * @param dt The time step for the simulation.
 */
void UpdateVelocity(std::vector<Particle>& particle, double dt) {
    for (auto& particleVelocity : particle) {
        // Acceleration = Force / Mass
        // ax = acceleration of x
        // ay = acceleration of y
        // az = acceleration of z
        double ax = particleVelocity.fx / particleVelocity.mass;
        double ay = particleVelocity.fy / particleVelocity.mass;
        double az = particleVelocity.fz / particleVelocity.mass;
        

        //velocity = velocity + force / mass * dt
        //acceleration = force / mass
        //velocity = velocity + acceleration * dt
        particleVelocity.vx = particleVelocity.vx + ax * dt;
        particleVelocity.vy = particleVelocity.vy + ay * dt;
        particleVelocity.vz = particleVelocity.vz + az * dt;

    }

}

/**
 * Updates the positions of all particles in the simulation.
 * This function is used to update the positions of all particles in the simulation
 * using the updated velocities from the previous step.
 * 
 * @param particles A vector of Particle objects representing all particles in the simulation.
 * @param dt The time step for the simulation.
 */
void UpdatePosition(std::vector<Particle>& particle, double dt) {
    for (auto& par : particle) {
        // Update the position 
        //position = position + velocity * dt
        par.x = par.x + par.vx * dt;
        par.y = par.y + par.vy * dt;
        par.z = par.z + par.vz * dt;
    }
}


/**
 * Outputs the state of the simulation to a file.
 * This function is used to output the state of the simulation to a file.
 * The state includes the mass, position, velocity, and force for each particle.
 * The file format is as follows:
 * 1. The number of particles in the simulation.
 * 2. For each particle: the mass, x position, y position, z position, x velocity, y velocity, z velocity, x force, y force, and z force.
 * 3. Each value is separated by a tab character. "\t"
 * 
 * @param particles A vector of Particle objects representing all particles in the simulation.
 * @param file The ofstream object representing the file to output the state to.
 */
void OutputState(const std::vector<Particle>& particles, std::ofstream& file) {
    file << particles.size() << std::endl;
    for (const auto& p : particles) {
        file << p.mass << "\t" << p.x << "\t" << p.y << "\t" << p.z << "\t" 
        << p.vx << "\t" << p.vy << "\t" << p.vz << "\t" 
        << p.fx << "\t" << p.fy << "\t" << p.fz << std::endl;
    }
}


int main(int argc, char* argv[]) {
    if(argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <num_particles> <dt> <num_steps> <dump_interval>" << std::endl;
        return 1;
    }
        int num_particles = std::stoi(argv[1]);
        double dt = std::stod(argv[2]);
        int num_steps = std::stoi(argv[3]);
        double dump_interval = std::stod(argv[4]);

    std::vector<Particle> particles(num_particles);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1000.0, 1000.0);
    for (auto& particle : particles) {
        particle.mass = 1.0;
        particle.x = dis(gen);
        particle.y = dis(gen);
        particle.z = dis(gen);
        particle.vx = 0.0;
        particle.vy = 0.0;
        particle.vz = 0.0;
    }

    std::ofstream file("solar.tsv");
    for (int i = 0; i < num_steps; i++) {
        CalculateForce(particles);
        UpdateVelocity(particles, dt);
        UpdatePosition(particles, dt);
        if (std::fmod(i, dump_interval) == 0) {
            OutputState(particles, file);
        }
    }
    file.close();
    return 0;
    /***
    std::ofstream file("small.tsv");
    for (int i = 0; i < num_steps; i++) {
        CalculateForce(particles);
        UpdateVelocity(particles, dt);
        UpdatePosition(particles, dt);
        if (i % dump_interval == 0) {
            OutputState(particles, file);
        }
    }
    file.close();
    return 0;
    ***/
}
