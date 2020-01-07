/*──────────────────────────────────────────────────────────────────────────── *
 *                                                                             *
 * |> DistNBody-HM ~ [Main] <|                                                 *
 *                                                                             *
 * MPI/OpenMP-parallel gravitational collisional N-body simulator via modern   *
 * C++ and hypermodularity.                                                    *
 *                                                                             *
 * (Maybe) updated version and other (related) cool stuff:                     *
 * https://github.com/emaballarin/distnbody-hm                                 *
 *                                                                             *
 *                                                                             *
 * Copyright (C) 2020-* Emanuele Ballarin <emanuele@ballarin.cc>.              *
 * All rights reserved. Distribution: Apache License 2.0.                      *
 * Full license text: https://ballarin.cc/legal/licenses/apacheii.txt          *
 *                                                                             *
 *                                                                             *
 * "Any fool can write code that a computer can understand. Good programmers   *
 *  write code that humans can understand."                                    *
 *                                               (M. Fowlers)                  *
 *                                                                             *
 * ─────────────────────────────────────────────────────────────────────────── */


/******************
 * GLOBAL PRAGMAS *
 ******************/

#pragma once


/*******************
 * GLOBAL INCLUDES *
 *******************/

#include <algorithm>  // STL algorithms on iterables (arrays, vectors)
#include <mpi.h>      // MPI parallelism
//#include <omp.h>        // OpenMP parallelism
#include <fcntl.h>      // File manipulations (POSIX)
#include <string>       // STL strings for filenames
#include <sys/stat.h>   // File manipulations (POSIX)
#include <sys/types.h>  // File manipulations (POSIX)
#include <unistd.h>     // File manipulations (POSIX)
#include <vector>       // STL Vectors


/******************
 * LOCAL INCLUDES *
 ******************/

#include "immutables.hpp"   // Mathematical constants
#include "physics.hpp"      // Physical evolution
#include "precision.hpp"    // Handling of real numbers
#include "triparticle.hpp"  // Objects that represent material bodies
#include "trivector.hpp"    // Objects that represent 3D vectors
#include "utilityfx.hpp"    // Utility functions


/***********
 * DEFINES *
 ***********/
#define _GNU_SOURCE                 // Explicitly enable GNU extensions
#define ull unsigned long long int  // A shorthand
#define SEED 35791246               // Seed for number generator
#define li long int
#define lli long long int
#define ui unsigned int
#define ul unsigned long int


/******************************************************************************
 ******************************************************************************/

/*********************
 * UTILITY FUNCTIONS *
 *********************/

inline real_t drand48_vel_adjust(real_t uniform0_1)
{
    return (uniform0_1 * 0.10) - 0.05;
}


/********
 * MAIN *
 ********/

int main(int argc, char* argv[])
{

    /* Define global, problem-specific constants */
    constexpr real_t GlobalMass = 100;
    constexpr ull Nparticles = 100;  // Change this if required!
    constexpr bool write_a_checkpoint{true};

    /* Compute global, problem-specific constants */
    constexpr real_t singlemass = GlobalMass / Nparticles;
    constexpr int Niter = 2;

    /* Define required MPI variables */
    int procId;
    int errCode;
    int P;
    MPI_Status status;
    int master{0};

    /* Initialize MPI */
    errCode = MPI_Init(&argc, &argv);
    errCode = MPI_Comm_size(MPI_COMM_WORLD, &P);
    errCode = MPI_Comm_rank(MPI_COMM_WORLD, &procId);

    /*
     * ~~ BRIEFING: ~~
     *
     * From now on:
     * - P is the nodes total;
     * - procId is the node index.
     *
     * This means that:
     * - Loops iterate in the range [0, P) or [0, P-1] --> USE THIS: for (auto i{0}; i < P; ++i)
     * - procId is the procId+1-th node.
     */

    /* GENERATE THE INITIAL CONDITIONS */

    // How many we need? (per node)
    ull particles_to_gen = workslice<ull>(Nparticles, static_cast<ull>(P), static_cast<ull>(procId)).stop
                           - workslice<ull>(Nparticles, static_cast<ull>(P), static_cast<ull>(procId)).start;

    // Seed the RNG
    srand48(SEED * (procId + 1));

    // Generate position/velocity of particles in arrays

    auto genXpos = new real_t[particles_to_gen];
    auto genXvel = new real_t[particles_to_gen];
    auto genYpos = new real_t[particles_to_gen];
    auto genYvel = new real_t[particles_to_gen];
    auto genZpos = new real_t[particles_to_gen];
    auto genZvel = new real_t[particles_to_gen];

    for (ull i = 0; i < particles_to_gen; ++i)
    {
        genXpos[i] = drand48();  // Already in the [0, 1) range!
        genYpos[i] = drand48();  // Already in the [0, 1) range!
        genZpos[i] = drand48();  // Already in the [0, 1) range!
        genXvel[i] = drand48_vel_adjust(drand48());
        genYvel[i] = drand48_vel_adjust(drand48());
        genZvel[i] = drand48_vel_adjust(drand48());
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /* Write the IC to files (one per node; no more complex IO is needed) */

    // Ordinary IO management
    std::string myFileName = std::to_string(procId);   // We will call files with the (unique) procId
    myFileName.append("_IC");                          // Better tag this, as the following will be different
    int myFileId = open(myFileName.c_str(),            // DO NOT USE NON-CONST METHODS ON myFileName FROM NOW ON!
                        O_WRONLY | O_CREAT | O_TRUNC,  // O_TRUNC -> O_APPEND if opening again
                        0644);                         // Change mode if needed by whatever reason

    /* Prepare & Dump data */

    // Precision of floats -> sizeof(<whatever has type real_t>)
    real_t writethis_realt{};
    int writethis = static_cast<int>(sizeof(writethis_realt));
    write(myFileId, &writethis, sizeof(writethis));

    // Particles in the file
    writethis = static_cast<int>(particles_to_gen);  // Change this when dumping during evolution
    write(myFileId, &writethis, sizeof(writethis));

    // The total number of files
    writethis = static_cast<int>(P);
    write(myFileId, &writethis, sizeof(writethis));

    // Time
    double writethis_double = 0.0;  // Change this when dumping during evolution
    write(myFileId, &writethis_double, sizeof(writethis_double));

    // Positions and velocities + Energy

    // Ooops, we need energies for that (kinetic only, fortunately!)
    auto Ekin_forIC = new real_t[particles_to_gen];

    for (ull particle_being_written{0}; particle_being_written < particles_to_gen; ++particle_being_written)
    {
        Ekin_forIC[particle_being_written] = singlemass
                                             * (genXvel[particle_being_written] * genXvel[particle_being_written]
                                                + genYvel[particle_being_written] * genYvel[particle_being_written]
                                                + genZvel[particle_being_written] * genZvel[particle_being_written])
                                             / 2.0;
    }


    // Rework this when dumping during evolution
    for (ull particle_being_written{0}; particle_being_written < particles_to_gen; ++particle_being_written)
    {
        // X
        writethis_realt = genXpos[particle_being_written];
        write(myFileId, &writethis_realt, sizeof(writethis_double));

        // Y
        writethis_realt = genYpos[particle_being_written];
        write(myFileId, &writethis_realt, sizeof(writethis_double));

        // Z
        writethis_realt = genYpos[particle_being_written];
        write(myFileId, &writethis_realt, sizeof(writethis_double));

        // vX
        writethis_realt = genXvel[particle_being_written];
        write(myFileId, &writethis_realt, sizeof(writethis_double));

        // vY
        writethis_realt = genYvel[particle_being_written];
        write(myFileId, &writethis_realt, sizeof(writethis_double));
        // vX
        writethis_realt = genZvel[particle_being_written];
        write(myFileId, &writethis_realt, sizeof(writethis_double));

        // Energies (kin)
        writethis_realt = Ekin_forIC[particle_being_written];
        write(myFileId, &writethis_realt, sizeof(writethis_double));
    }


    // Close what we opened
    close(myFileId);

    // Clean up
    delete[] Ekin_forIC;
    ////////////////////////////////////////////////////////////////////////////////////////////////


    /* DECOMPOSE THE DOMAIN (in a smart way!) */

    /*
     * In the following section we will broadcast ALL POSITIONS AND VELOCITIES from/to every node.
     * This may seem inefficient (or downright stupid), and therefore we need to write here...
     *
     * ~~ AN EXPLANATION: ~~
     *
     * Suppose that you have already sent ALL POSITIONS via a broadcast to/from every node already
     * (they will be needed for dynamical evolution anyway).
     *
     * Now you need to:
     *
     * 1) In each node, compute how many particles VELOCITIES to receive from each other node;
     * 2) Send/receive such velocities;
     * 3) Organize everything in such way that is useful for the dynamical evolution computation.
     *
     * In order to fulfil (2) we can choose to:
     * a) Broadcast velocities as we did for positions; then prune velocities in receiving node. This is what we will do.
     *    -> This means broadcasting N elements with a single broadcast at O(log(P))
     *       -> Complexity of: ~O(N*log(P))
     *
     * b) Prune velocities in sending node; then send/receive with 2-staggered parallel Send/Receive (optimal solution).
     *    -> This means sending ~(N/P) elements
     *                                          |> P/2 times if P is even
     *                                          |> (3/2)*(P-1) times if P is odd
     *       -> Complexity (avg.) of: O(N)
     *
     * Obviously, the (b) solution has to be preferred if the ideal case is considered and efficiency alone has
     * to be taken into account. This idea, also, assumes constant bandwidth among all nodes and P >> 4.
     * It also imples a less robust implementation (which has to be hand-coded by us).
     *
     * For that reasons, we deem adequate (though not optimal) to pursue path (1). This way, at the cost of
     * a slight (initial: this operation happens only once!) overhead allows to be fully topology-agnostic
     * (Squyres, 2007) and less bug-prone.
     */

    // We just need to know the amount of particles that given node has generated

    ull generated_by_node[P];  // Requires C >= 99

    for (int i = 0; i < P; ++i)
    {
        if (i == procId)
        {
            generated_by_node[i] = particles_to_gen;
        }
        else
        {
            generated_by_node[i] = workslice<ull>(Nparticles, static_cast<ull>(P), static_cast<ull>(i)).stop
                                   - workslice<ull>(Nparticles, static_cast<ull>(P), static_cast<ull>(i)).start;
        }
    }


    // We now need to receive the positions/velocities from all other nodes.
    // We need also to store them in (a) suitable structure(s), to give order to everything!

    /*
     * ~~ A trick, right-off the 'tricks hat'... ~~
     */
    std::vector<std::vector<real_t>> genPosMatX(P);
    std::vector<std::vector<real_t>> genPosMatY(P);
    std::vector<std::vector<real_t>> genPosMatZ(P);
    std::vector<std::vector<real_t>> genVelMatX(P);
    std::vector<std::vector<real_t>> genVelMatY(P);
    std::vector<std::vector<real_t>> genVelMatZ(P);

    for (int i = 0; i < P; i++)
    {
        genPosMatX[i].resize(generated_by_node[i]);
        genPosMatY[i].resize(generated_by_node[i]);
        genPosMatZ[i].resize(generated_by_node[i]);
        genVelMatX[i].resize(generated_by_node[i]);
        genVelMatY[i].resize(generated_by_node[i]);
        genVelMatZ[i].resize(generated_by_node[i]);
    }
    /*
     * In a different context, this should be considered to be idiotic. However, in this case...
     *
     * Here we need to build a "thing" that has:
     * - 2 dimensions;
     * - Both dimensions known only at runtime;
     * - Preferrably one dimension has to dependent on the other (i.e. a list of non-equidimensional arrays)
     * - Preferrably each row has to be castable to and from raw arrays.
     *
     * std::vectors are slower than arrays (not too much, though!) and this seems the most reasonable
     * alternative to boost::multiarrays, std::lists of boost::arrays of other similar results
     * of type-mumbo-jumbo.
     *
     * The alternative is to write a custom class in order to that. It is surely the best thing to do.
     * However, it will take time, a lot of testing and a sparkle of genial attitude to build
     * something more battle-tested than C++'s STL library.
     *
     * Probably fair enough.
     */

    /* BROADCAST POSITIONS */

    for (int whichNode{0}; whichNode < P; ++whichNode)
    {
        auto mpibuffer = new real_t[generated_by_node[whichNode]];

        if (whichNode == procId)
        {
            /* Store in `generatedPos/VelMatrix[whichNode]` and broadcast to all */

            // Xpos
            mpibuffer = genXpos;  // Copy to avoid overwriting
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genPosMatX[whichNode].begin());
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);

            // Ypos
            mpibuffer = genYpos;  // Copy to avoid overwriting
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genPosMatY[whichNode].begin());
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);

            // Zpos
            mpibuffer = genZpos;  // Copy to avoid overwriting
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genPosMatZ[whichNode].begin());
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);

            // Xvel
            mpibuffer = genXvel;  // Copy to avoid overwriting
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genVelMatX[whichNode].begin());
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);

            // Yvel
            mpibuffer = genYvel;  // Copy to avoid overwriting
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genVelMatY[whichNode].begin());
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);

            // Zvel
            mpibuffer = genZvel;  // Copy to avoid overwriting
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genVelMatZ[whichNode].begin());
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
        }
        else
        {
            /* Receive from node `whichNode` and store in `generatedPos/VelMatrix[whichNode]` */

            // Xpos
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genPosMatX[whichNode].begin());

            // Ypos
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genPosMatY[whichNode].begin());

            // Zpos
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genPosMatZ[whichNode].begin());

            // Xvel
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genVelMatX[whichNode].begin());

            // Yvel
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genVelMatY[whichNode].begin());

            // Zvel
            MPI_Bcast(mpibuffer, generated_by_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
            std::copy(mpibuffer, mpibuffer + generated_by_node[whichNode], genVelMatZ[whichNode].begin());
        }

        // Be easy on memory...
        delete[] mpibuffer;
    }

    // Now we have every stored the generated position/velocities (X, Y, Z) on every node...
    // ... and we can free some memory (NOTE: `gen[_]pos/vel` is now filled with garbage anyway)!
    delete[] genXpos;
    delete[] genYpos;
    delete[] genZpos;
    delete[] genXvel;
    delete[] genYvel;
    //delete[] genZvel;   // Has already been deleted by `delete[] mpibuffer` and some STL magic ;)

    /*
     * Now we need to:
     * - Prune unneeded velocities;
     * - Order everything in structures more suitable to dynamical evolution and
     *   communication of updated positions.
     */

    // Local storage of positions (per-thread) from now on
    std::vector<std::vector<real_t>> LocalPositionsX(P);
    std::vector<std::vector<real_t>> LocalPositionsY(P);
    std::vector<std::vector<real_t>> LocalPositionsZ(P);

    // Local storage of velocities from now on
    // (Don't abuse with velocities memory...)
    std::vector<real_t> LocalVelocitiesX(0);
    std::vector<real_t> LocalVelocitiesY(0);
    std::vector<real_t> LocalVelocitiesZ(0);

    // Don't abuse with position memory either...
    for (int i = 0; i < P; i++)
    {
        LocalPositionsX[i].resize(0);
        LocalPositionsY[i].resize(0);
        LocalPositionsZ[i].resize(0);
    }

    /* Perform some looping... ;) */

    for (int whichNode = 0; whichNode < P; ++whichNode)  // Over nodes
    {
        for (ull whichParticle = 0; whichParticle < generated_by_node[whichNode]; ++whichParticle)  // Over particles
        {
            // Compute assignee node
            real_t currentX = genPosMatX[whichNode].at(whichParticle);
            int assignee = static_cast<int>(((currentX) / (1.0 / P)));  // 1.0 is the box side length

            // Place in the right place of the new matrix
            LocalPositionsX[assignee].push_back(currentX);
            LocalPositionsY[assignee].push_back(genPosMatY[whichNode].at(whichParticle));
            LocalPositionsZ[assignee].push_back(genPosMatZ[whichNode].at(whichParticle));

            // Prune velocities and add to matrix if needed
            if (assignee == procId)
            {
                LocalVelocitiesX.push_back(genVelMatX[whichNode].at(whichParticle));
                LocalVelocitiesY.push_back(genVelMatY[whichNode].at(whichParticle));
                LocalVelocitiesZ.push_back(genVelMatZ[whichNode].at(whichParticle));
            }
        }
    }

    // Now we have what we wanted and unused std::vectors are automatically deleted through smart pointers! ;)

    // For the sake of practicality...
    ull particles_of_node[P];                            // Requires C >= 99
    for (int whichNode = 0; whichNode < P; ++whichNode)  // Over nodes
    {
        particles_of_node[whichNode] = LocalPositionsX[whichNode].size();
    }

    /*
     * ~~ BRIEFING: ~~
     *
     * At this point we find ourselves with a structure that can handle the dynamical evolution of the system
     * and at the same time allow for an easy synchronization of positions among threads. Before being able to
     * implement the "dynamical evolution" loop we just need to load positions and velocities of a given assigned
     * node in a purposefully-built TriParticle object, which will remain local to each node.
     *
     * Dynamical evolution of each particle will happen thanks to built-in methods of a TriParticle; synchronization
     * will be realized by broadcasting the positions vector composed of the X, Y, Z coordinates of the TriParticles
     * after the evolution is computed.
     *
     * This (seemingly) unnecessary data copy allows a build-and-swap style of thread safety, as the broadcast can only
     * take place after all the evolutions have been computed.
     *
     * Due to the fact that the expected number of particles-per-node is the same among nodes, we tolerate minimal
     * work unbalance among nodes, and we choose not to address specifically the load-balancing problem.
     */

    // In each node, only before first iteration...
    auto LocalTriParticles = new TriPart[particles_of_node[procId]];

    for (ull thisParticle = 0; thisParticle < particles_of_node[procId]; ++thisParticle)
    {
        // For each particle this node needs to evolve...
        LocalTriParticles[thisParticle].m(singlemass);
        LocalTriParticles[thisParticle].x(LocalPositionsX[procId].at(thisParticle));
        LocalTriParticles[thisParticle].y(LocalPositionsY[procId].at(thisParticle));
        LocalTriParticles[thisParticle].z(LocalPositionsZ[procId].at(thisParticle));
        LocalTriParticles[thisParticle].vx(LocalVelocitiesX.at(thisParticle));
        LocalTriParticles[thisParticle].vy(LocalVelocitiesY.at(thisParticle));
        LocalTriParticles[thisParticle].vz(LocalVelocitiesZ.at(thisParticle));
        LocalTriParticles[thisParticle].E(LocalTriParticles[thisParticle].Ekin());
    }

    // Now we have what we wanted and unused std::vectors are automatically deleted through smart pointers! ;)
    // e.g. velocities not in TriParticles.


    real_t simulation_time = 0.0;
    for (int iteration{0}; iteration <= Niter + 1; ++iteration)
    {
        /* Now we fully unroll the first dynamical iteration... */

        // REMEMBER TO FREE MEMORY!

        // This will be required to keep track of the forces acting on particles
        auto LocalForces_acting_on = new TriVec[particles_of_node[procId]];

        // This will be required to keep track of the potential energies of PREVIOUS step (more on this later...)
        auto LocalPotentialEnergy_of = new real_t[particles_of_node[procId]];

        // This will be required to compute the new timestep
        real_t maxmodulus{0.0};

        // For each particle assigned to the node...
        for (ull thisParticle{0}; thisParticle < particles_of_node[procId]; ++thisParticle)
        {
            TriVec ForceOn_thisParticle{0, 0, 0};  // The to-be force acting on such particle
            real_t EnergyOf_thisParticle = 0.0;

            // For each particle in the system (from the local position matrix)...
            for (int whichNode{0}; whichNode < P; ++whichNode)
            {
                for (ull whichParticle{0}; whichParticle < particles_of_node[whichNode]; ++whichParticle)
                {
                    TriVec deltaForce = {0, 0, 0};
                    real_t deltaEnergy = 0;
                    if (whichNode != procId || thisParticle != whichParticle)
                    {
                        // If particles are different, compute mutual interaction...
                        deltaForce = Fgrav_on1_from2(
                          singlemass, singlemass,
                          {LocalPositionsX[procId].at(thisParticle), LocalPositionsY[procId].at(thisParticle),
                           LocalPositionsZ[procId].at(thisParticle)},
                          {LocalPositionsX[whichNode].at(whichParticle), LocalPositionsY[whichNode].at(whichParticle),
                           LocalPositionsZ[whichNode].at(whichParticle)});

                        deltaEnergy = Egrav_of1_from2(
                          singlemass, singlemass,
                          {LocalPositionsX[procId].at(thisParticle), LocalPositionsY[procId].at(thisParticle),
                           LocalPositionsZ[procId].at(thisParticle)},
                          {LocalPositionsX[whichNode].at(whichParticle), LocalPositionsY[whichNode].at(whichParticle),
                           LocalPositionsZ[whichNode].at(whichParticle)});
                    }
                    ForceOn_thisParticle += deltaForce;
                    EnergyOf_thisParticle += deltaEnergy;
                }
            }

            // Keep track of computed forces (which are invariant)
            LocalForces_acting_on[thisParticle] = ForceOn_thisParticle;
            LocalPotentialEnergy_of[thisParticle] = EnergyOf_thisParticle;


            // Compute candidate max modulus
            real_t candidate_maxmodulus =
              singlemass * (mod(ForceOn_thisParticle) / mod(LocalTriParticles[thisParticle].vel()));

            // Update global max modulus if needed
            if (candidate_maxmodulus > maxmodulus)
            {
                maxmodulus = candidate_maxmodulus;
            }
        }

        // Compute ENERGIES OF PREVIOUS ITERATION!
        for (ull thisParticle = 0; thisParticle < particles_of_node[procId]; ++thisParticle)
        {
            LocalTriParticles[thisParticle].E(LocalTriParticles[thisParticle].Ekin()
                                              + LocalPotentialEnergy_of[thisParticle]);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////
        // HERE WE HAVE ALL THE DATA WE NEED FOR PREVIOUS ITERATION! //

        if (write_a_checkpoint && iteration != 0)  // We already wrote the conditions for iteration = -1 (i.e. the IC)
        {  // Also: NEVER simplfy the expression: write_a_checkpoint CAN BE constexpr!
            ///

            // Ordinary IO management
            myFileName = std::to_string(procId);           // We will call files with the (unique) procId
            myFileName.append("_ECP");                     // ECP = Evolution CheckPoint
            myFileId = open(myFileName.c_str(),            // DO NOT USE NON-CONST METHODS ON myFileName FROM NOW ON!
                            O_WRONLY | O_CREAT | O_TRUNC,  // Overwrite previous (but not if _IC)
                            0644);                         // Change mode if needed by whatever reason

            /* Prepare & Dump data */

            // Precision of floats -> sizeof(<whatever has type real_t>)
            writethis_realt = 0;
            writethis = static_cast<int>(sizeof(writethis_realt));
            write(myFileId, &writethis, sizeof(writethis));

            // Particles in the file
            writethis = static_cast<int>(particles_of_node[procId]);
            write(myFileId, &writethis, sizeof(writethis));

            // The total number of files
            writethis = static_cast<int>(P);
            write(myFileId, &writethis, sizeof(writethis));

            // Time
            writethis_double = simulation_time;
            write(myFileId, &writethis_double, sizeof(writethis_double));

            // Positions and velocities + Energy
            for (ull particle_being_written{0}; particle_being_written < particles_of_node[procId]; ++particle_being_written)
            {
                // X
                writethis_realt = LocalTriParticles[particle_being_written].x();
                write(myFileId, &writethis_realt, sizeof(writethis_double));

                // Y
                writethis_realt = LocalTriParticles[particle_being_written].y();
                write(myFileId, &writethis_realt, sizeof(writethis_double));

                // Z
                writethis_realt = LocalTriParticles[particle_being_written].z();
                write(myFileId, &writethis_realt, sizeof(writethis_double));

                // vX
                writethis_realt = LocalTriParticles[particle_being_written].vx();
                write(myFileId, &writethis_realt, sizeof(writethis_double));

                // vY
                writethis_realt = LocalTriParticles[particle_being_written].vy();
                write(myFileId, &writethis_realt, sizeof(writethis_double));
                // vX
                writethis_realt = LocalTriParticles[particle_being_written].vz();
                write(myFileId, &writethis_realt, sizeof(writethis_double));

                // Energies (kin)
                writethis_realt = LocalTriParticles[particle_being_written].E();
                write(myFileId, &writethis_realt, sizeof(writethis_double));
            }


            // Close what we opened
            close(myFileId);
        }

        // IF EXITING
        if (iteration == Niter+1)
        {
            delete[] LocalForces_acting_on;
            delete[] LocalPotentialEnergy_of;
            MPI_Finalize();
            return 0;
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////


        // Compute the global max modulus
        real_t global_maxmodulus{0.0};
        MPI_Allreduce(&maxmodulus, &global_maxmodulus, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Compute new timestep (max timestep)
        real_t new_timestep = veleps / global_maxmodulus;

        // Execute the dynamic step on all particles...
        for (ull thisParticle = 0; thisParticle < particles_of_node[procId]; ++thisParticle)
        {
            LocalTriParticles[thisParticle].stepForce_CBBC(LocalForces_acting_on[thisParticle], new_timestep,
                                                           1.0);  // Box BC
        }

        // Free some memory
        delete[] LocalForces_acting_on;
        delete[] LocalPotentialEnergy_of;


        /*
         * WHAT WE HAVE NOW?
         * - The `LocalTriParticles` vector contains updated positions and updated velocities
         * - The `Local[Positions|Velocities][X, Y, Z]` need to be updated:
         *                                              |> Positions both locally and remotely
         *                                              |> Velocities only locally
         */


        // Prepare the position X, Y, Z vectors to be sent/updated locally

        auto updatedXpos = new real_t[particles_of_node[procId]];
        auto updatedYpos = new real_t[particles_of_node[procId]];
        auto updatedZpos = new real_t[particles_of_node[procId]];

        for (ull particle = 0; particle < particles_of_node[procId]; ++particle)
        {
            updatedXpos[particle] = LocalTriParticles[particle].x();
            updatedYpos[particle] = LocalTriParticles[particle].y();
            updatedZpos[particle] = LocalTriParticles[particle].z();
        }


        // (Update & Send) what needs to be sent
        for (int whichNode{0}; whichNode < P; ++whichNode)
        {
            auto mpibuffer = new real_t[particles_of_node[whichNode]];

            if (whichNode == procId)
            {
                // Xpos
                mpibuffer = updatedXpos;  // Copy to avoid overwriting
                std::copy(mpibuffer, mpibuffer + particles_of_node[whichNode], LocalPositionsX[whichNode].begin());
                MPI_Bcast(mpibuffer, particles_of_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);

                // Ypos
                mpibuffer = updatedYpos;  // Copy to avoid overwriting
                std::copy(mpibuffer, mpibuffer + particles_of_node[whichNode], LocalPositionsY[whichNode].begin());
                MPI_Bcast(mpibuffer, particles_of_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);

                // Zpos
                mpibuffer = updatedZpos;  // Copy to avoid overwriting
                std::copy(mpibuffer, mpibuffer + particles_of_node[whichNode], LocalPositionsZ[whichNode].begin());
                MPI_Bcast(mpibuffer, particles_of_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
            }
            else
            {
                // Xpos
                MPI_Bcast(mpibuffer, particles_of_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
                std::copy(mpibuffer, mpibuffer + particles_of_node[whichNode], LocalPositionsX[whichNode].begin());

                // Ypos
                MPI_Bcast(mpibuffer, particles_of_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
                std::copy(mpibuffer, mpibuffer + particles_of_node[whichNode], LocalPositionsY[whichNode].begin());

                // Zpos
                MPI_Bcast(mpibuffer, particles_of_node[whichNode], MPI_DOUBLE, whichNode, MPI_COMM_WORLD);
                std::copy(mpibuffer, mpibuffer + particles_of_node[whichNode], LocalPositionsZ[whichNode].begin());
            }
            // Be easy on memory...
            delete[] mpibuffer;
        }

        // Free some memory
        delete[] updatedXpos;
        delete[] updatedYpos;
        // delete[] updatedZpos;    // Has already been deleted by `delete[] mpibuffer` and some STL magic ;)

        // Update time (which will only be needed at next iteration for dumping to file)
        simulation_time += new_timestep;
    }


    // The following two lines should never execute, if the program is used correctly!
    MPI_Finalize();
    return 0;
}
