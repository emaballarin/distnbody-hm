# DistNBody-HM
An *MPI-parallel* (optionally hybrid: in-node *OpenMP-parallel*) implementation
of a *gravitational*, *collisional* N-body simulator, developed in *Modern C++*
according to the principles of *hypermodular development*.  

Please note that **the simulator here implemented is not physically accurate**
and that this is not a priority for such project.  

Parallelism and extreme modularity are the main goal. In addition to that,
however, only minimal effort is required in order to extend the simulator to
additional cases and specific requirements (at least in the case of
*conservative potentials*).

---
