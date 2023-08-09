// July 30th, 2023.
//
// Source code for the replica exchange (parallel tempering) method used in
// "Quantum-inspired Encoding Enhances Monte Carlo Sampling of Soft Matter Systems" by S. Longo, P. Hauke, P. Faccioli and C. Micheletti
//
// The code is provided  "as is" and with CC-BY-NC-SA license, with no liability nor implied assistance in compiling the code, running it etc.
//
// The code uses the ran2 random number generator of "Numerical Recipes" by Press et al., which must be provided by the user- or substitited with an equivalent routine/
//
// The code was developed by CM.
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_memory.c"
#include <string.h>
#include "ran2.c"
#include "data_structures.c"


#define EPS (1.0e-6) // threshold to detect if parameters are different from zero
#define N_replicas (8) // Number of Parallel Tempering replicas

/*  global variables */
int Lattice_Dimensionality; // 2 for square lattice, 3 for cubic
int L[3] ;  // Number of *sites* per each lattice edge, obviously we need *at most* 3
int Nbonds; // Number of bonds in polymer chain
int Ncorners; // Number of corners in polymer chain

long int idum= (-1); // Seed for NR random number gnerator ran2

#include "lattice_definitions.c" /* Contains construct_lattice_structure,
check_if_sites_are_neighbours and deallocate_memory XXX */
#include "mc_routines.c"
#include "io.c" // IO operations



int main(){

  struct Lattice lattice;
  struct System_properties system[N_replicas], trial_system;
  int i, j, k, l, m, dice;
  int Output_trj, Output_energy, Output_spin_config;
  long int Swap_interval, Dump_interval, N_dumps;
  long int n_moves, Ntot_moves, n_accepted;

  printf("Replica exchange (parallel tempering) sampling. Number of Replicas: %d\n\n",N_replicas);
  read_MC_parameters("INPUT_PT_PARAMETERS.DAT", &Swap_interval, &Dump_interval, &N_dumps, &Output_trj, &Output_energy, &Output_spin_config);
  read_input_lattice_parameters("INPUT_LATTICE_PARAMS.DAT");
  read_input_temperatures("INPUT_PT_REPLICA_TEMPERATURES.DAT",system);

//  printf("\tAllow_4loops: %3d\n",Allow_4loops);

  // Prepare the lattice
  construct_lattice_structure(&lattice);

  for(i=0; i < N_replicas; i++){
    system[i].n_spins=lattice.n_bonds;
    system[i].spin_array=c1t(system[i].n_spins);
    initialization(&lattice, &(system[i]), Output_trj, Output_energy, Output_spin_config);
  }
  trial_system.spin_array=c1t(system[0].n_spins);

  n_accepted=0;
  Ntot_moves=Dump_interval*N_dumps;
  for(n_moves=1; n_moves <= Ntot_moves; n_moves++){

      /**** NORMAL MOVES ****/
    for(i=0; i < N_replicas; i++){
      MC_plaquette_flip_move(&lattice, &(system[i]), &trial_system);
    }

      /**** SWAP MOVES ****/
    if((n_moves% Swap_interval==0)){
      dice = (int) floor(ran2(&idum)*2); /* start randomly from replica 0 or 1, */
      for(j=dice; j < (N_replicas-1); j+=2){
        if (j>=(N_replicas-1)) continue; // It should not happen, but better safe than sorry
        system[j].n_accepted_swaps+=attempt_swap(&system[j],&system[j+1], &trial_system);
        system[j].n_swaps++;
      }
    }


      if (n_moves%Dump_interval==0){
         for(i=0; i < N_replicas; i++){
          dump_observables(&lattice, &(system[i]), Output_trj, Output_energy, Output_spin_config);
        }
      }

      if (n_moves%(Ntot_moves/10)==0){
       for(i=0; i < N_replicas; i++){
        printf("\t Replica %3d  Accepted moves: %9ld out of %9ld  Accepted swaps: %8ld out of %8ld  current_energy: %4d\n",i,system[i].n_accepted_moves, n_moves, system[i].n_accepted_swaps, system[i].n_swaps, system[i].energy);
        }
        printf("\n");
      }

    }

// Exit gracefully.
  free_c1t(trial_system.spin_array);
  deallocate_memory(&lattice, system);
}


