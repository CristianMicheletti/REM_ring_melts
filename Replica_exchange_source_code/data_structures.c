/* Lattice is a struct variable that work as a hash table to helps calling the correct
variables in the main script (e.g. the indices of a site in lattice or in the qubit array,
the sites that flanks a bond, the bonds that make a triplet and so on */
struct Lattice {
  int n_sites;
  int **site;   // first index is for site enumeration, the second is for the Cartesian coordinates.
  int n_bonds;
  int **bond; // first index is for bond enumeration, the second is for the indices of the two sites.
  int n_corners;
  int **corner_bond; // first index is for triplet enumeration, the second is for the indices of the two bonds.
  int n_plaquettes;
  int **plaquette_bond; // first index is for plaquette enumeration, the second is for the indices of the four bonds, (indices 0-1 and 2-3 are for the two pairs of opposite bonds
  int **plaquette_corner; // first index is for plaquette enumeration, the second is for the indices of the four corners
};

struct System_properties{
  int n_spins;
  char *spin_array;
  int energy; // the energy is equal to the number of corners, which is proportional to the total curvature
  double temperature;
  long int n_accepted_moves;
  long int n_accepted_swaps;
  long int n_swaps;
  char trj_filename[100];
  char energy_filename[100];
  char configs_filename[100];
};
