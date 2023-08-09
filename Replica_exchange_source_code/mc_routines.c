int count_active_corners_of_plaquette(struct Lattice *lattice,  struct System_properties *system, int plaquette){

  int i, j, k, corner, n_active_corners;

  n_active_corners=0;

  for (i=0; i < 4; i++){
    corner=lattice->plaquette_corner[plaquette][i];
    // get the indices of the bonds forming corner i of the plaquette
    j = lattice->corner_bond[corner][0];
    k = lattice->corner_bond[corner][1];
    // corner is active if both bonds are active (=1)
    n_active_corners+=system->spin_array[j]*system->spin_array[k];
  }

  return(n_active_corners);
}


/****************************/

int count_all_active_corners(struct Lattice *lattice,  struct System_properties *system){

  int i, n_active_corners;
  int temp;

  n_active_corners=0;

  for (i=0; i< lattice->n_plaquettes; i++){
    // get the indices of the bonds forming corner i of the plaquette
  //    printf("Plaquette %d has %d corners\n",i,count_active_corners_of_plaquette(lattice,system,i));
    n_active_corners+=count_active_corners_of_plaquette(lattice,system,i);
  }
  return(n_active_corners);
}


/****************************/

int check_configuration_is_OK(struct Lattice *lattice,  struct System_properties *system){

  int i, j, site, c;
  char *site_occupation;

    site_occupation=c1t(lattice->n_sites);
    for(i=0; i < lattice->n_sites; i++){
      site_occupation[i]=0;
    }

    c=0;
    for(i=0; i < system->n_spins; i++){
//      printf("%d ",system->spin_array[i]);
      c+=system->spin_array[i];
    }
//    printf("\n");
//    printf("Number of active bonds: %d\n",c);

    for(i=0; i < system->n_spins; i++){
      if (system->spin_array[i]==1){
        for(j=0; j < 2; j++){
        site=lattice->bond[i][j];
        site_occupation[site]++;
        if (site_occupation[site]>2) {
          free_c1t(site_occupation);
//          printf("degree > 2\n");
          return(0);
        }
        }
      }
    }

// We have already excluded the cases of sites with three incident bonds
// Now let us exclude cases with only one incident bond

    for(i=0; i < lattice->n_sites; i++){
      if (site_occupation[i]==1) {
//        printf("degree = 1\n");
        free_c1t(site_occupation);
        return(0);
      }
    }
    free_c1t(site_occupation);
    return(1);
}

/***************************/
void initialization(struct Lattice *lattice,  struct System_properties *system, int Output_trj, int Output_energy, int Output_spin_config){

  int c, i, j;
  char filename[100];
  FILE *fp;

  printf("Initializing system.\n");


    // check if there exists a file containing an initial configuration for the considered lattice size.

  sprintf(filename,"INPUT_INITIAL_CONFIG.DAT");
  if ((fp = fopen (filename, "r")) == NULL) {
    printf("WARNING: Could not find input file %s.\n",filename);
    goto generate_config_stochastically;
  }
    // File exists, let's use it!

  // first check that the file has the correct number of entries

  c=0;
  while(fscanf(fp,"%d",&i)!=EOF){
    if ((i!=0) && (i!=1)){
      printf("FATAL ERROR. Entry in %s is different from 0/1  (%d)\n",filename,i);
      exit(1);
    }
    c++;
  }
  if (c!=system->n_spins){
    printf("FATAL ERROR. File %s has wrong number of entries (%d vs %d)\n",filename,c,system->n_spins);
    exit(1);
  }

  rewind(fp);
  for(i=0; i < system->n_spins; i++){
    fscanf(fp,"%d",&j);
    system->spin_array[i]=j;
  }
  printf("Read initial configuration from %s\n",filename);
  goto initialize_variables;


 generate_config_stochastically:
  printf("Generating initial configuration stochastically\n");
  do{
    for(i=0; i < system->n_spins; i++) system->spin_array[i]=0;



    // otherwise generate the initial configuration randomly, this works only for very small systems
    c=0;
    do{
      j= (int) floor(ran2(&idum)*lattice->n_bonds);
      if (system->spin_array[j]==0){
        system->spin_array[j]=1;
        c++;
      }
    }while (c < Nbonds);

    }while (check_configuration_is_OK(lattice,system)==0);

 initialize_variables:
  system->energy=count_all_active_corners(lattice,system);
  system->n_accepted_moves=0;
  system->n_accepted_swaps=0;
  system->n_swaps=0;

  sprintf(system->trj_filename,"output_trj_T%e.XYZ",system->temperature);
  sprintf(system->energy_filename,"output_energy_T%e.dat",system->temperature);
  sprintf(system->configs_filename,"output_spin_configs_T%e.dat",system->temperature);

  if (Output_trj==1){
    fp=fopen(system->trj_filename,"w");
    fclose(fp);
  }
  if (Output_energy==1){
    fp=fopen(system->energy_filename,"w");
    fclose(fp);
  }
  if (Output_spin_config==1){
    fp=fopen(system->configs_filename,"w");
    fclose(fp);
  }

  printf("Completed initialization at full filling.\n");
  printf("Initial energy: %d\n",system->energy);
}

/***************************/

void initialization_from_scratch(struct Lattice *lattice,  struct System_properties *system){
 /* practical only for very small systems */

  /*   initialize_system_at_full_filling(system,Lattice); */

  int c, i, j;
  FILE *fp;

  printf("Initializing system.\n");
  do{
    for(i=0; i < system->n_spins; i++) system->spin_array[i]=0;


    // random configuration
    c=0;
    do{
      j= (int) floor(ran2(&idum)*lattice->n_bonds);
      if (system->spin_array[j]==0){
        system->spin_array[j]=1;
        c++;
      }
    }while (c < Nbonds);

    }while (check_configuration_is_OK(lattice,system)==0);

  system->energy=count_all_active_corners(lattice,system);
  system->n_accepted_moves=0;
  system->n_accepted_swaps=0;
  system->n_swaps=0;

  sprintf(system->trj_filename,"trj_T%e.XYZ",system->temperature);
  sprintf(system->energy_filename,"energy_T%e.dat",system->temperature);

  fp=fopen(system->trj_filename,"w");
  fclose(fp);
  fp=fopen(system->energy_filename,"w");
  fclose(fp);

  printf("Completed initialization at full filling.\n");
}

/******************************************/

void copy_system_a_to_b(struct System_properties *system_a, struct System_properties *system_b){

  int i;


  system_b->n_spins=system_a->n_spins;
  for(i=0; i < system_a->n_spins; i++){
    system_b->spin_array[i]=system_a->spin_array[i];
  }
  system_b->energy=system_a->energy;
  system_b->temperature=system_a->temperature;
  system_b->n_accepted_moves=system_a->n_accepted_moves;
  system_b->n_accepted_swaps=system_a->n_accepted_swaps;
  system_b->n_swaps=system_a->n_swaps;
}


/***************************/

void MC_plaquette_flip_move(struct Lattice *lattice,  struct System_properties *system, struct System_properties *trial_system){

  int i, j, k, bond[4], n_flips, DeltaE;
  int flip_flag;


  copy_system_a_to_b(system,trial_system);



  if (ran2(&idum)<0.5){
    n_flips=1; // 50/50 probability of a single flip, local change
  }
  else{
    n_flips= 2+ (int) floor(ran2(&idum)*9); // 50/50 probability of doing multiple flips, from 2 to 10.
  }

 flip_flag=1; // initialization of the flip flag */

  for(k=0; k < n_flips; k++){
    /* pick random plaquette */
    i = (int) floor(ran2(&idum)*lattice->n_plaquettes);

  /* check if plaquette is flippable */

    for(j=0; j < 4; j++){
      bond[j]=lattice->plaquette_bond[i][j];
    }

  /* Plaquette is flippable in two cases only:
      (a) spin_bond0=spin_bond1=1 and spin_bond2=spin_bond3=0
      (b) spin_bond0=spin_bond1=0 and spin_bond2=spin_bond3=1
    */

  // case (a)
    if ( ((trial_system->spin_array[bond[0]] + trial_system->spin_array[bond[1]]) ==2) &&
      ((trial_system->spin_array[bond[2]] + trial_system->spin_array[bond[3]]) ==0) ){
    /* flip it */
      trial_system->spin_array[bond[0]]=0;
      trial_system->spin_array[bond[1]]=0;
      trial_system->spin_array[bond[2]]=1;
      trial_system->spin_array[bond[3]]=1;
      trial_system->energy= count_all_active_corners(lattice, trial_system);
    }

 // case (b)
    else if ( ((trial_system->spin_array[bond[0]] + trial_system->spin_array[bond[1]]) ==0) &&
      ((trial_system->spin_array[bond[2]] + trial_system->spin_array[bond[3]]) ==2) ){
    /* flip it */
      trial_system->spin_array[bond[0]]=1;
      trial_system->spin_array[bond[1]]=1;
      trial_system->spin_array[bond[2]]=0;
      trial_system->spin_array[bond[3]]=0;
      trial_system->energy=count_all_active_corners(lattice, trial_system);
    }
    else flip_flag=0; // plaquette was not flippable
  }

  // if flip_flag=0, then the series of attempted flips was not successful.
  //  Reject the move, i.e. return without doing anything.

  if (flip_flag==1){
    // series of attempted flips were ok, let's invoke the Metropolis criterion
    DeltaE=system->energy-trial_system->energy;
    if (log(ran2(&idum))< (DeltaE/system->temperature)){
      // We passed the Metropolis criterion, accept the move
       copy_system_a_to_b(trial_system,system);
       system->n_accepted_moves++;
    }
  }
}


/*******************/

int attempt_swap(struct System_properties  *system1, struct System_properties  *system2,  struct System_properties *trial_system){

  int j, temp, swapped;
  long int litemp;
  double old_weight, new_weight, temp2;


  old_weight=system1->energy/system1->temperature+system2->energy/system2->temperature;
  new_weight=system1->energy/system2->temperature+system2->energy/system1->temperature;

  swapped=0;
  if ( (old_weight - new_weight ) > log(ran2(&idum))){
    swapped=1;
    /* Let's swap the two systems */

    copy_system_a_to_b(system1,trial_system);
    copy_system_a_to_b(system2,system1);
    copy_system_a_to_b(trial_system,system2);

    /* but fix back the temperatures and info on number of moves*/
    temp2 =system1->temperature;
    system1->temperature=system2->temperature;
    system2->temperature=temp2;

    litemp =system1->n_accepted_moves;
    system1->n_accepted_moves=system2->n_accepted_moves;
    system2->n_accepted_moves=litemp;

    litemp =system1->n_accepted_swaps;
    system1->n_accepted_swaps=system2->n_accepted_swaps;
    system2->n_accepted_swaps=litemp;

    litemp =system1->n_swaps;
    system1->n_swaps=system2->n_swaps;
    system2->n_swaps=litemp;

  }

  return(swapped);
}
