void construct_lattice_structure(struct Lattice *lattice);
int check_if_sites_are_neighbours(int i, int j, struct Lattice *lattice);
int check_if_bonds_are_in_same_plaquette(int i, int j, struct Lattice *lattice);
int manhattan_distance(int i, int j, struct Lattice *lattice);
double d_distance(int i, int j, struct Lattice *lattice);
int check_if_bonds_are_opposite(int i, int j, struct Lattice *lattice);
int check_if_sites_are_a_straight_or_corner_triplet(int i, int j, int k, struct Lattice *lattice, int *corner_flag);
void deallocate_memory(struct Lattice *lattice, struct System_properties *system);
int find_corner_index_for_two_incident_bonds(struct Lattice *lattice, int bond1, int bond2);

/*******************************************************************************/


void construct_lattice_structure(struct Lattice *lattice){

  int i, j, k, l, m, coord[Lattice_Dimensionality], counter, corner_flag;
  int bond1, bond2;



  /* Establish coordinates of sites of D-dimensional lattice. */

  printf("Mapping sites, bonds and corners of %d-dimensional lattice of sides ",Lattice_Dimensionality);
  for(i=0; i < Lattice_Dimensionality; i++) {printf(" %d",L[i]);}
  printf("\n");

  lattice->n_sites=1;
  for(i=0; i < Lattice_Dimensionality; i++){
    lattice->n_sites*=L[i];
  }
  lattice->site= i2t(lattice->n_sites,Lattice_Dimensionality);

  for(j=0; j < Lattice_Dimensionality; j++){
    coord[j]=0;
  }

  for(i=0; i < lattice->n_sites; i++){
    for(j=0; j < Lattice_Dimensionality; j++){
      lattice->site[i][j]= coord[j];
    }

    k=0;
    coord[k]++;
    while (coord[k]==L[k]) {
      coord[k]=0;
      k++;
      coord[k]++;
    };
  }

  printf("Constructed Cartesian array of %d lattice sites.\n",lattice->n_sites);
  Nbonds=lattice->n_sites;
  printf("Set Nbonds = %d - full filling condition required for applying bond flips\n", Nbonds);
/*  We need to pre-enumerate bond and corners, in order to allocate just enough memory space.
    For lattices with simple geometries this pre-enumeration could be done analytically, but we do it
    explicitly so the code can be adapted to e.g. different boundary conditions or particular constraints. */


/* Pre-enumerating and then constructing hash table array of distinct (oriented) bonds */
  counter=0;
  for(i=0; i < lattice->n_sites; i++){
    for(j=i+1; j < lattice->n_sites; j++){
       if (check_if_sites_are_neighbours(i,j, lattice)==1) counter++;
    }
  }

  lattice->n_bonds=counter;
  lattice->bond= i2t(lattice->n_bonds,2);
  counter=0;
  for(i=0; i < lattice->n_sites; i++){
    for(j=i+1; j < lattice->n_sites; j++){
       if (check_if_sites_are_neighbours(i,j, lattice)==1) {
          lattice->bond[counter][0]=i;
          lattice->bond[counter][1]=j;
          counter++;
       }
    }
  }

if (counter != lattice->n_bonds){
    printf("Fatal error, inconsistent counting of bonds: %d vs %d\n",counter,lattice->n_bonds);
    exit(0);
  }

  printf("Constructed hash table of %d distinct oriented bonds.\n",lattice->n_bonds);



/* Pre-enumerating and then constructing hash table array of distinct (oriented) triplets, i, j, k,
  with j being the middle site, and i <k . Notice that index j is NOT restricted.

  We enumerate corner triplets only.

 */


/* check corner triplets */

  counter=0;
  for(l=0; l < lattice->n_bonds; l++){
    for(m=l+1; m < lattice->n_bonds; m++){
      // bond must be distinct but must share one site, j, the middle site of the triplet.
      if (lattice->bond[m][0]==lattice->bond[l][0] ) {j=lattice->bond[l][0]; i=lattice->bond[l][1]; k=lattice->bond[m][1];}
      else if (lattice->bond[m][1]==lattice->bond[l][0] ) {j=lattice->bond[l][0]; i=lattice->bond[l][1]; k=lattice->bond[m][0];}
      else if (lattice->bond[m][0]==lattice->bond[l][1] ) {j=lattice->bond[l][1]; i=lattice->bond[l][0]; k=lattice->bond[m][1];}
      else if (lattice->bond[m][1]==lattice->bond[l][1] ) {j=lattice->bond[l][1]; i=lattice->bond[l][0]; k=lattice->bond[m][0];}
      else continue;
      if (check_if_sites_are_a_straight_or_corner_triplet(i,j,k,lattice, &corner_flag)==1) {
        /* i, j, k is a triplet. Now check if it is a straight or corner one and select as needed */
          if (corner_flag==1) counter++;
      }
    }
  }


  lattice->n_corners=counter;
  lattice->corner_bond= i2t(lattice->n_corners,2);

  counter=0;
  for(l=0; l < lattice->n_bonds; l++){
    for(m=l+1; m < lattice->n_bonds; m++){
      // bond must be distinct but must share one site, j, the middle site of the triplet.
      if (lattice->bond[m][0]==lattice->bond[l][0] ) {j=lattice->bond[l][0]; i=lattice->bond[l][1]; k=lattice->bond[m][1];}
      else if (lattice->bond[m][1]==lattice->bond[l][0] ) {j=lattice->bond[l][0]; i=lattice->bond[l][1]; k=lattice->bond[m][0];}
      else if (lattice->bond[m][0]==lattice->bond[l][1] ) {j=lattice->bond[l][1]; i=lattice->bond[l][0]; k=lattice->bond[m][1];}
      else if (lattice->bond[m][1]==lattice->bond[l][1] ) {j=lattice->bond[l][1]; i=lattice->bond[l][0]; k=lattice->bond[m][0];}
      else continue;

      if (check_if_sites_are_a_straight_or_corner_triplet(i,j,k,lattice, &corner_flag)==1) {
        /* i, j, k is a triplet. Now check if it is a straight or corner one and select as needed */
        if (corner_flag==1) {
          lattice->corner_bond[counter][0]=l;
          lattice->corner_bond[counter][1]=m;
          counter++;
        }
      }

    }
  }

 if (counter != lattice->n_corners){
    printf("Fatal error, inconsistent counting of corners: %d vs %d\n",counter,lattice->n_corners);
    exit(0);
  }
  printf("Constructed hash table of %d distinct corner triplets\n",lattice->n_corners);



  /* Same scheme to pre-enumerate plaquettes */

  counter=0;
  for(i=0; i < lattice->n_bonds; i++){
    for(j=i+1; j < lattice->n_bonds; j++){
      if (check_if_bonds_are_in_same_plaquette(i,j, lattice)==0) continue;
      for(k=j+1; k < lattice->n_bonds; k++){
        if (check_if_bonds_are_in_same_plaquette(i,k, lattice)==0) continue;
        if (check_if_bonds_are_in_same_plaquette(j,k, lattice)==0) continue;
        for(l=k+1; l < lattice->n_bonds; l++){
          if (check_if_bonds_are_in_same_plaquette(i,l, lattice)==0) continue;
          if (check_if_bonds_are_in_same_plaquette(j,l, lattice)==0) continue;
          if (check_if_bonds_are_in_same_plaquette(k,l, lattice)==0) continue;
          counter++;
        }
      }
    }
  }

  lattice->n_plaquettes=counter;

  lattice->plaquette_bond=i2t(lattice->n_plaquettes,4);
  counter=0;
  for(i=0; i < lattice->n_bonds; i++){
    for(j=i+1; j < lattice->n_bonds; j++){
      if (check_if_bonds_are_in_same_plaquette(i,j, lattice)==0) continue;
      for(k=j+1; k < lattice->n_bonds; k++){
        if (check_if_bonds_are_in_same_plaquette(i,k, lattice)==0) continue;
        if (check_if_bonds_are_in_same_plaquette(j,k, lattice)==0) continue;
        for(l=k+1; l < lattice->n_bonds; l++){
          if (check_if_bonds_are_in_same_plaquette(i,l, lattice)==0) continue;
          if (check_if_bonds_are_in_same_plaquette(j,l, lattice)==0) continue;
          if (check_if_bonds_are_in_same_plaquette(k,l, lattice)==0) continue;

          // three options for opposite bonds in the plaquette:
          // (a) i-j and k-l ;
          // (b) i-k and j-l ;
          // (c) i-l and j-k .

            if (check_if_bonds_are_opposite(i,j,lattice)==1) {
              // option a
              lattice->plaquette_bond[counter][0]=i;
              lattice->plaquette_bond[counter][1]=j;
              lattice->plaquette_bond[counter][2]=k;
              lattice->plaquette_bond[counter][3]=l;
              counter++;
            }


            else if (check_if_bonds_are_opposite(i,k,lattice)==1) {
              // option b
              lattice->plaquette_bond[counter][0]=i;
              lattice->plaquette_bond[counter][1]=k;
              lattice->plaquette_bond[counter][2]=j;
              lattice->plaquette_bond[counter][3]=l;
              counter++;
            }


            else if (check_if_bonds_are_opposite(i,l,lattice)==1) {
              // option l
              lattice->plaquette_bond[counter][0]=i;
              lattice->plaquette_bond[counter][1]=l;
              lattice->plaquette_bond[counter][2]=j;
              lattice->plaquette_bond[counter][3]=k;
              counter++;
            }
        }
      }
    }
  }

  if (counter != lattice->n_plaquettes){
    printf("Fatal error, inconsistent counting of plaquettes: %d vs %d\n",counter,lattice->n_plaquettes);
    exit(1);

  }
  printf("Constructed hash table of bonds %d distinct plaquettes.\n",lattice->n_plaquettes);


  lattice->plaquette_corner=i2t(lattice->n_plaquettes,4);


  for(i=0; i < lattice->n_plaquettes; i++){
    // we have four pairs of bonds invoved in corners: (0) 0-2,;(1) 0-3; (2) 1-2; (3) 1-3

    for (k=0; k < 4 ; k++){

      switch(k){

      case 0:
        bond1=lattice->plaquette_bond[i][0];
        bond2=lattice->plaquette_bond[i][2];
        break;

      case 1:
        bond1=lattice->plaquette_bond[i][0];
        bond2=lattice->plaquette_bond[i][3];
        break;

      case 2:
        bond1=lattice->plaquette_bond[i][1];
        bond2=lattice->plaquette_bond[i][2];
        break;

      case 3:
        bond1=lattice->plaquette_bond[i][1];
        bond2=lattice->plaquette_bond[i][3];
        break;
      }

    find_corner_index_for_two_incident_bonds(lattice, bond1, bond2);
    lattice->plaquette_corner[i][k]=find_corner_index_for_two_incident_bonds(lattice, bond1, bond2);
  }

  }

  printf("Constructed hash table of corners of %d distinct plaquettes.\n",lattice->n_plaquettes);


}

/********************/

int manhattan_distance(int i, int j, struct Lattice *lattice){

  int l, delta, manhattan_dist_ij;

  manhattan_dist_ij=0;
  for(l=0; l < Lattice_Dimensionality; l++){
    delta=abs(lattice->site[i][l] - lattice->site[j][l]);
    manhattan_dist_ij+=delta;
  }

  return(manhattan_dist_ij);
}

/********************/

double d_distance(int i, int j, struct Lattice *lattice){

  int l;
  double delta;
  double dist_ij;

  dist_ij=0;
  for(l=0; l < Lattice_Dimensionality; l++){
    delta=abs(lattice->site[i][l] - lattice->site[j][l]);
    dist_ij+=delta*delta;
  }
  dist_ij=sqrt(dist_ij);

  return(dist_ij);
}

/********************/

int check_if_sites_are_neighbours(int i, int j, struct Lattice *lattice){


  if (manhattan_distance(i,j,lattice) != 1) return(0);
  return(1);
}

/********************/


int check_if_bonds_are_in_same_plaquette(int i, int j, struct Lattice *lattice){

  int site1_bond1, site2_bond1,  site1_bond2, site2_bond2;

    site1_bond1=lattice->bond[i][0];
    site2_bond1=lattice->bond[i][1];

    site1_bond2=lattice->bond[j][0];
    site2_bond2=lattice->bond[j][1];

    if (d_distance(site1_bond1,site1_bond2,lattice) > (sqrt(2) + EPS)) return(0);
    else if (d_distance(site1_bond1,site2_bond2,lattice) > (sqrt(2) + EPS)) return(0);
    else if (d_distance(site2_bond1,site1_bond2,lattice) > (sqrt(2) + EPS)) return(0);
    else if (d_distance(site2_bond1,site2_bond2,lattice) > (sqrt(2) + EPS)) return(0);
    else return (1);

}

/********************/

int check_if_bonds_are_opposite(int i, int j, struct Lattice *lattice){

   int site1_bond1, site2_bond1,  site1_bond2, site2_bond2;

    site1_bond1=lattice->bond[i][0];
    site2_bond1=lattice->bond[i][1];

    site1_bond2=lattice->bond[j][0];
    site2_bond2=lattice->bond[j][1];

    if (d_distance(site1_bond1,site1_bond2,lattice) < EPS) return(0);
    else if (d_distance(site1_bond1,site2_bond2,lattice) < EPS) return(0);
    else if (d_distance(site2_bond1,site1_bond2,lattice) < EPS) return(0);
    else if (d_distance(site2_bond1,site2_bond2,lattice) < EPS) return(0);
    else return(1);
}

/*********************/


int check_if_sites_are_a_straight_or_corner_triplet(int i, int j, int k, struct Lattice *lattice, int *corner_flag){

 int l, delta, manhattan_dist_ij, manhattan_dist_jk, manhattan_dist_ik ;
 double Eucl_dist_ik;

 /*
 printf("** Testing %d %d %d\n",i,j,k);
 l=i;
 printf("Coord of site %d   %2d %2d\n",l,lattice->site[l][0],lattice->site[l][1]);
 l=j;
 printf("Coord of site %d   %2d %2d\n",l,lattice->site[l][0],lattice->site[l][1]);
 l=k;
 printf("Coord of site %d   %2d %2d\n",l,lattice->site[l][0],lattice->site[l][1]);
 */

 /* i and j must be at Manhattan distance 1 */
  manhattan_dist_ij=0;
  for(l=0; l < Lattice_Dimensionality; l++){
    delta=abs(lattice->site[i][l] - lattice->site[j][l]);
    manhattan_dist_ij+=delta;
  }
  //  printf("Manhattan distance_ij: %d\n",manhattan_dist_ij);
  if (manhattan_dist_ij != 1) return(0);

/* j and k must be at Manhattan distance 1 */
  manhattan_dist_jk=0;
  for(l=0; l < Lattice_Dimensionality; l++){
    delta=abs(lattice->site[j][l] - lattice->site[k][l]);
    manhattan_dist_jk+=delta;
  }
  //  printf("Manhattan distance_jk: %d\n",manhattan_dist_jk);
  if (manhattan_dist_jk != 1) return(0);

/* i and k must be at Manhattan distance 2 */
  manhattan_dist_ik=0;
  Eucl_dist_ik=0.0;
  for(l=0; l < Lattice_Dimensionality; l++){
    delta=abs(lattice->site[i][l] - lattice->site[k][l]);
    manhattan_dist_ik+=delta;
    Eucl_dist_ik+=delta*delta;
  }
  //  printf("Manhattan distance_ik: %d  sq_Eucl_distance: %lf\n",manhattan_dist_ik,Eucl_dist_ik);
  if (manhattan_dist_ik != 2) return(0);
  Eucl_dist_ik=sqrt(Eucl_dist_ik);

/* i,j,k are a viable triplet. Now check if it is a corner or straight triplet */

  *corner_flag=0;
  if (Eucl_dist_ik < (manhattan_dist_ik*0.99999999)) *corner_flag=1;

  //  printf("Corner flag: %d\n",*corner_flag);
  return(1);

}

/********************/

int find_corner_index_for_two_incident_bonds(struct Lattice *lattice, int bond1, int bond2){

 int i, j, k;

 for(i=0; i <lattice->n_corners; i++){
  j=lattice->corner_bond[i][0];
  k=lattice->corner_bond[i][1];
  if (((bond1==j) && (bond2==k)) || ((bond1==k) && (bond2==j)) ){
    return(i);
  }
 }

 printf("Fatal Error. Incident bonds not found in list of corners\n");
 exit(1);

}

/********************/
void deallocate_memory(struct Lattice *lattice, struct System_properties *system){

  int i;
  printf("De-allocating memory.\n");


  for(i=0; i < N_replicas; i++){
  free_c1t(system[i].spin_array);
  }
  free_i2t(lattice->plaquette_corner);
  free_i2t(lattice->plaquette_bond);
  free_i2t(lattice->bond);
  free_i2t(lattice->site);
}
