/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [TriParticle] <|                                          *
 *                                                                             *
 * An efficient and user-friendly structure for 3D material points             *
 * (a.k.a. particles)                                                          *
 *                                                                             *
 * (Maybe) updated version and other (related) cool stuff:                     *
 * https://github.com/emaballarin/distnbody-hm                                 *
 *                                                                             *
 *                                                                             *
 * Copyright (C) 2020-* Emanuele Ballarin <emanuele@ballarin.cc>.              *
 * All rights reserved. Distribution: Apache License 2.0.                      *
 * Full license text: https://ballarin.cc/legal/licenses/apacheii.txt          *
 *                                                                             *
 * ─────────────────────────────────────────────────────────────────────────── */


/******************
 * GLOBAL PRAGMAS *
 ******************/

#pragma once


/************
 * INCLUDES *
 ************/

#include <cassert>
#include "precision.hpp"  // Handling of real numbers
#include "utilityfx.hpp"  // Utility functions
#include "trivector.hpp"  // 3D vectors
#include "physics.hpp"    // Physical evolution

/******************************************************************************
 ******************************************************************************/


class TriPart
{
    public:

    /********
    * CTORS *
    *********/

    inline TriPart() : position{}, velocity{}, energy{}, mass{}
    {
    }

    inline explicit TriPart(real_t _m) : position{}, velocity{}, energy{}, mass(_m)
    {
        assert(_m > 0); // Mass must be positive
    }

    inline TriPart(TriVec _pos, real_t _m) : position(_pos), velocity{}, energy{}, mass(_m)
    {
        assert(_m > 0);  // Mass must be positive
    }

    inline TriPart(TriVec _pos, TriVec _vel, real_t _m) : position(_pos), velocity(_vel), energy{}, mass(_m)
    {
        assert(_m > 0);  // Mass must be positive
    }

    // NOTE: Eneegy initialization MUST be done manually after allocation! (otw: E=0)

    /***********
     * GETTERS *
     ***********/

    inline TriVec pos() const
    {
        return position;
    }

    inline real_t x() const
    {
        return position.x;
    }

    inline real_t y() const
    {
        return position.y;
    }

    inline real_t z() const
    {
        return position.z;
    }

    inline TriVec vel() const
    {
        return velocity;
    }

    inline real_t vx() const
    {
        return velocity.x;
    }

    inline real_t vy() const
    {
        return velocity.y;
    }

    inline real_t vz() const
    {
        return velocity.z;
    }

    inline real_t E() const
    {
        return energy;
    }

    inline real_t m() const
    {
        return mass;
    }


    /***********
     * SETTERS *
     ***********/

    inline void pos(TriVec _pos)
    {
        position.x = _pos.x;
        position.y = _pos.y;
        position.z = _pos.z;
    }

    inline void pos(real_t _pos[3]) // Do not const-ify!
    {
        /* Do not assign null pointer dereferences */
        assert(_pos != nullptr);
        assert(_pos + 1 != nullptr);
        assert(_pos + 2 != nullptr);

        position.x = _pos[0];
        position.y = _pos[1];
        position.z = _pos[2];
    }

    inline void x(real_t _x)
    {
        position.x = _x;
    }

    inline void y(real_t _y)
    {
        position.y = _y;
    }

    inline void z(real_t _z)
    {
        position.z = _z;
    }

    inline void vel(TriVec _vel)
    {
        position.x = _vel.x;
        position.y = _vel.y;
        position.z = _vel.z;
    }

    inline void vel(real_t _vel[3])  // Do not const-ify!
    {
        /* Do not assign null pointer dereferences */
        assert(_vel != nullptr);
        assert(_vel + 1 != nullptr);
        assert(_vel + 2 != nullptr);

        position.x = _vel[0];
        position.y = _vel[1];
        position.z = _vel[2];
    }

    inline void vx(real_t _vx)
    {
        velocity.x = _vx;
    }

    inline void vy(real_t _vy)
    {
        velocity.y = _vy;
    }

    inline void vz(real_t _vz)
    {
        velocity.z = _vz;
    }

    inline void E(real_t _E)
    {
        assert(_E >= 0);
        energy = _E;
    }

    inline void m(real_t _mass)
    {
        assert(_mass >= 0);
        mass = _mass;
    }

   /**************
    * MEMBER FXS *
    **************/

    inline real_t Ekin() const
   {
        return (mass*dot(velocity, velocity))/2;
   }

   // Impose box boundary conditions for position, in-place
   inline void boxBC(real_t boxx, real_t boxy, real_t boxz)
   {
       assert(boxx >= 0);
       assert(boxy >= 0);
       assert(boxz >= 0);

       position.x = rangeloop<long int>(position.x, boxx);
       position.y = rangeloop<long int>(position.y, boxy);
       position.z = rangeloop<long int>(position.z, boxz);
   }

   // Impose CUBIC box boundary conditions for position, in-place
   inline void boxBC_cubic(real_t boxside)
   {
       assert(boxside >= 0);

       position.x = rangeloop<long int>(position.x, boxside);
       position.y = rangeloop<long int>(position.y, boxside);
       position.z = rangeloop<long int>(position.z, boxside);
   }

   /* Physical evolutors */

   inline void stepForce(TriVec _force, real_t _timestep)
   {
       TriVec _newvel = NewVel(velocity, accel_from_Fm(_force, mass), _timestep);
       velocity.x = _newvel.x;
       velocity.y = _newvel.y;
       velocity.z = _newvel.z;

       TriVec _newpos = NewPos(position, _newvel, _timestep);
       position.x = _newpos.x;
       position.y = _newpos.y;
       position.z = _newpos.z;
   }

   inline void stepForce_CBBC(TriVec _force, real_t _timestep, real_t _boxside)
   {
       this->stepForce(_force, _timestep);
       this->boxBC_cubic(_boxside);
   }

    private:
    TriVec position;
    TriVec velocity;
    real_t energy;
    real_t mass;
    //const real_t mass;
};
