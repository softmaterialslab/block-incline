#ifndef _PARTICLE_H
#define _PARTICLE_H

#include"vector3d.h"

class PARTICLE{
    //members
public:
    double m;           // mass of the particle
    int ty;             // type of the particle
    VECTOR3D position;  // position of the particle in 3D
    VECTOR3D velocity;  // velocity of the particle in 3D
    VECTOR3D force;     // force on the particle in 3D
    double ke;          // kinetic energy
    double pe;          // potential energy


    //make a particle
    PARTICLE(double initial_m = 1, int initial_ty = 1, VECTOR3D initial_position = VECTOR3D(0,0,0), VECTOR3D initial_velocity = VECTOR3D(0,0,0))
    {
        m = initial_m;
        ty = initial_ty;
        position = initial_position;
        velocity = initial_velocity;
    }

    // member function
    void kinetic_energy()
    {
        ke = 0.5 * m * velocity.Magnitude() * velocity.Magnitude();
    }

    // code member function to update position

    void update_position(double dt)
    {
        position = position + (velocity * dt);
    }

    // code member function to update velocity

    void update_velocity(double half_dt)
    {
        velocity = velocity + (force * (half_dt/m));
    }

};

#endif
