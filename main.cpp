
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
// 1. Physical constants
constexpr double G = 6.67430e-11;
constexpr double AU = 1.496e11; // Astronomical Unit in meters

// 2. Simulation settings
constexpr double DT = 60.0;
constexpr int NUM_STEPS =  5256000;

// 3. Data structures
struct Body {
    std::string name;
    double mass;
    double position[3];
    double velocity[3];
    double acceleration[3];
    double potential_energy;
    double kinetic_energy;
};
struct Derivative {
    double velocity[3];
    double acceleration[3];
};

// 4. Functions
double calculate_potential_energy(const std::vector<Body>& bodies) {
    double potential_energy = 0.0;

    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i + 1; j < bodies.size(); ++j) {

            double dx = bodies[i].position[0] - bodies[j].position[0];
            double dy = bodies[i].position[1] - bodies[j].position[1];
            double dz = bodies[i].position[2] - bodies[j].position[2];

            double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (distance == 0.0) continue;

            potential_energy -= G * bodies[i].mass * bodies[j].mass / distance;
        }
    }

    return potential_energy;
}
void calculate_acceleration(std::vector<Body>& bodies) {
    // Gravitational acceleration calculations will go here
    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].acceleration[0] = 0.0;
        bodies[i].acceleration[1] = 0.0;
        bodies[i].acceleration[2] = 0.0;
        for (size_t j = 0; j < bodies.size(); ++j) {
            if (i != j) {
                double dx = bodies[j].position[0] - bodies[i].position[0];
                double dy = bodies[j].position[1] - bodies[i].position[1];
                double dz = bodies[j].position[2] - bodies[i].position[2];
                double distance_squared = dx * dx + dy * dy + dz * dz;
                if (distance_squared == 0) continue;
                double distance = std::sqrt(distance_squared);
                double acceleration_magnitude = G * bodies[j].mass / distance_squared;
                bodies[i].acceleration[0] += acceleration_magnitude * (dx / distance) ;
                bodies[i].acceleration[1] += acceleration_magnitude * (dy / distance) ;
                bodies[i].acceleration[2] += acceleration_magnitude * (dz / distance) ;
            }
        }
    }
}   
void semi_euler_update(std::vector<Body>& bodies) {
    for (auto& body : bodies) {
        body.velocity[0] += body.acceleration[0] * DT;
        body.velocity[1] += body.acceleration[1] * DT;
        body.velocity[2] += body.acceleration[2] * DT;

        body.position[0] += body.velocity[0] * DT;
        body.position[1] += body.velocity[1] * DT;
        body.position[2] += body.velocity[2] * DT;
    }
}
void euler_update(std::vector<Body>& bodies) {
    for (auto& body : bodies) {
        body.position[0] += body.velocity[0] * DT;
        body.position[1] += body.velocity[1] * DT;
        body.position[2] += body.velocity[2] * DT;

        body.velocity[0] += body.acceleration[0] * DT;
        body.velocity[1] += body.acceleration[1] * DT;
        body.velocity[2] += body.acceleration[2] * DT;

        body.kinetic_energy = 0.5 * body.mass * (body.velocity[0] * body.velocity[0] + body.velocity[1] * body.velocity[1] + body.velocity[2] * body.velocity[2]);
        
    }
    double potential_energy = calculate_potential_energy(bodies);

    for (auto& body : bodies) {
        body.potential_energy = potential_energy;
    }
}
std::vector<Derivative> evaluate(const std::vector<Body>& bodies,
                                 const std::vector<Derivative>& d,
                                 double dt);

void RK4_updater(std::vector<Body>& bodies){
    auto original = bodies;
    std::vector<Derivative> initial_derivatives;
    calculate_acceleration(original);
    for (size_t i = 0; i < original.size(); ++i) {
        initial_derivatives.push_back({
            {0, 0, 0},
            {0, 0, 0}
        });
    }
    auto k1 = evaluate(original, initial_derivatives, DT);
    auto k2 = evaluate(original, k1, DT / 2);
    auto k3 = evaluate(original, k2, DT / 2);
    auto k4 = evaluate(original, k3, DT);

    for (size_t i = 0; i < bodies.size(); ++i) {
        bodies[i].position[0] = original[i].position[0] +
            (k1[i].velocity[0] + 2.0 * k2[i].velocity[0] +
             2.0 * k3[i].velocity[0] + k4[i].velocity[0]) * DT / 6.0;
        bodies[i].position[1] = original[i].position[1] +
            (k1[i].velocity[1] + 2.0 * k2[i].velocity[1] +
             2.0 * k3[i].velocity[1] + k4[i].velocity[1]) * DT / 6.0;
        bodies[i].position[2] = original[i].position[2] +
            (k1[i].velocity[2] + 2.0 * k2[i].velocity[2] +
             2.0 * k3[i].velocity[2] + k4[i].velocity[2]) * DT / 6.0;

        bodies[i].velocity[0] = original[i].velocity[0] +
            (k1[i].acceleration[0] + 2.0 * k2[i].acceleration[0] +
             2.0 * k3[i].acceleration[0] + k4[i].acceleration[0]) * DT / 6.0;
        bodies[i].velocity[1] = original[i].velocity[1] +
            (k1[i].acceleration[1] + 2.0 * k2[i].acceleration[1] +
             2.0 * k3[i].acceleration[1] + k4[i].acceleration[1]) * DT / 6.0;
        bodies[i].velocity[2] = original[i].velocity[2] +
            (k1[i].acceleration[2] + 2.0 * k2[i].acceleration[2] +
             2.0 * k3[i].acceleration[2] + k4[i].acceleration[2]) * DT / 6.0;

    }
    calculate_acceleration(bodies);
    double potential_energy = calculate_potential_energy(bodies);

    for (auto& body : bodies) {
        body.potential_energy = potential_energy;
    }
}

std::vector<Derivative> evaluate(const std::vector<Body>& bodies,const std::vector<Derivative>& d, double dt){
   std::vector<Body> temporary = bodies;
    for (size_t i = 0; i < temporary.size(); ++i) {
        temporary[i].position[0] += d[i].velocity[0] * dt;
        temporary[i].position[1] += d[i].velocity[1] * dt;
        temporary[i].position[2] += d[i].velocity[2] * dt;
        temporary[i].velocity[0] += d[i].acceleration[0] * dt;
        temporary[i].velocity[1] += d[i].acceleration[1] * dt;
        temporary[i].velocity[2] += d[i].acceleration[2] * dt;
    }
    calculate_acceleration(temporary);
    std::vector<Derivative> result;
    for (size_t i = 0; i < temporary.size(); ++i) {
        result.push_back({
    {temporary[i].velocity[0], temporary[i].velocity[1], temporary[i].velocity[2]},
    {temporary[i].acceleration[0], temporary[i].acceleration[1], temporary[i].acceleration[2]}});
    }
    return result;
}

bool write_history_to_file(const std::vector<std::vector<Body>>& body_history, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return false;
    }
    file << "step,body,mass,x,y,z,vx,vy,vz,ax,ay,az,kinetic_energy,potential_energy\n";
    for (size_t step = 0; step < body_history.size(); ++step) {
        for (const auto& body : body_history[step]) {
            file << step << "," 
                 << body.name << "," 
                 << body.mass << ","
                 << body.position[0] << "," 
                 << body.position[1] << "," 
                 << body.position[2] << ","
                 << body.velocity[0] << ","
                 << body.velocity[1] << ","
                 << body.velocity[2] << ","
                 << body.acceleration[0] << ","
                 << body.acceleration[1] << ","
                 << body.acceleration[2] << ","
                 << body.kinetic_energy << ","
                 << body.potential_energy << "\n";
        }
    }

    file.close();
    return true;
}
// 5. Main simulation
int main() {
     // Store history for one body (e.g., Earth)
    Body earth;
    earth.name = "Earth";
    earth.mass = 5.972e24;
    earth.position[0] = 1.496e11;
    earth.position[1] = 0.0;
    earth.position[2] = 0.0;
    earth.velocity[0] = 0.0;
    earth.velocity[1] = 29780.0;
    earth.velocity[2] = 0.0;
    earth.acceleration[0] = 0.0;
    earth.acceleration[1] = 0.0;
    earth.acceleration[2] = 0.0;
    Body sun;
    sun.name = "Sun";
    sun.mass = 1.989e30;
    sun.position[0] = 0.0;
    sun.position[1] = 0.0;
    sun.position[2] = 0.0;
    sun.velocity[0] = 0.0;
    sun.velocity[1] = 0.0;
    sun.velocity[2] = 0.0;
    sun.acceleration[0] = 0.0;
    sun.acceleration[1] = 0.0;
    sun.acceleration[2] = 0.0;
    double potential_energy = -G * earth.mass * sun.mass / sqrt((earth.position[0] - sun.position[0]) * (earth.position[0] - sun.position[0]) + (earth.position[1] - sun.position[1]) * (earth.position[1] - sun.position[1]) + (earth.position[2] - sun.position[2]) * (earth.position[2] - sun.position[2]));
    sun.kinetic_energy = 0.5 * sun.mass * (sun.velocity[0] * sun.velocity[0] + sun.velocity[1] * sun.velocity[1] + sun.velocity[2] * sun.velocity[2]);
    sun.potential_energy = potential_energy;
    earth.potential_energy = potential_energy;
    earth.kinetic_energy = 0.5 * earth.mass * (earth.velocity[0] * earth.velocity[0] + earth.velocity[1] * earth.velocity[1] + earth.velocity[2] * earth.velocity[2]);
    std::vector<Body> initial_bodies = {earth, sun};
    std::vector<Body> bodies = initial_bodies;
    std::vector<std::vector<Body>> body_history_euler;
    body_history_euler.push_back({bodies});
    for (int step = 0; step < NUM_STEPS; ++step) {
        calculate_acceleration(bodies);
        euler_update(bodies);

        body_history_euler.push_back({bodies});
    };
    write_history_to_file(body_history_euler, "body_history_euler.csv");
    bodies = initial_bodies;
    std::vector<std::vector<Body>> body_history_rk4;
    body_history_rk4.push_back({bodies});
    for (int step = 0; step < NUM_STEPS; ++step) {
   
        RK4_updater(bodies);

        body_history_rk4.push_back({bodies});
    };

    write_history_to_file(body_history_rk4, "body_history_rk4.csv");
}