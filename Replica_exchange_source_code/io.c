void read_MC_parameters(char *filename, long int *Swap_interval, long int *Dump_interval,long int *N_dumps, int *Output_trj, int *Output_energy, int *Output_spin_config){

  FILE *fp;
  char line[200], descriptor[100], value[100];

// Default parameters
 *Swap_interval=0;      // required parameter
 *Dump_interval=0;      // required parameter
 *N_dumps=0;            // required parameter
 *Output_trj=1;         // optional parameter
 *Output_energy=1;      // optional parameter
 *Output_spin_config=1; // optional parameter

 printf("Reading MC parameters from input file %s\n",filename);
 fp=fopen(filename,"r");

 while (fgets(line, 200, fp)){
    sscanf(line,"%s %s",descriptor,value);
    printf("%s %s\n",descriptor,value);

    if (strcmp (descriptor, "Swap_interval")==0){
      sscanf(value,"%ld",Swap_interval);
    }

  if (strcmp (descriptor, "Dump_interval")==0){
      sscanf(value,"%ld",Dump_interval);
    }

    if (strcmp (descriptor, "N_dumps")==0){
      sscanf(value,"%ld",N_dumps);
    }

    if (strcmp (descriptor, "Output_trj")==0){
      sscanf(value,"%d",Output_trj);
      if ((*Output_trj!=0) && (*Output_trj!=1)){
        printf("FATAL ERROR: %s can only be 0 or 1\n",descriptor);
        exit(1);
      }
    }

    if (strcmp (descriptor, "Output_energy")==0){
      sscanf(value,"%d",Output_energy);
      if ((*Output_energy!=0) && (*Output_energy!=1)){
        printf("FATAL ERROR: %s can only be 0 or 1\n",descriptor);
        exit(1);
      }
    }

    if (strcmp (descriptor, "Output_spin_config")==0){
      sscanf(value,"%d",Output_spin_config);
      if ((*Output_spin_config!=0) && (*Output_spin_config!=1)){
        printf("FATAL ERROR: %s can only be 0 or 1\n",descriptor);
        exit(1);
      }
    }
  }



  if (*Swap_interval<1){
    printf("FATAL ERROR: Swap_interval <1 or missing in file %s\n",filename);
    exit(1);
  }
  if (*Dump_interval<1){
    printf("FATAL ERROR: Dump_interval <1 or missing in file %s\n",filename);
    exit(1);
  }

  if (*N_dumps<1){
    printf("FATAL ERROR: N_dumps <1 or missing in file %s\n",filename);
    exit(1);
  }

  printf("Finished reading MC parameters\n\n");
}

/************************/

void read_input_lattice_parameters(char *filename){

  FILE *fp;
  char line[200], descriptor[100], value[100];
  int i, n_box_sides;

  n_box_sides=0;
  printf("Reading lattice parameters from input file %s\n",filename);
  fp=fopen(filename,"r");

  while (fgets(line, 200, fp)){
    sscanf(line,"%s %s",descriptor,value);
    printf("%s %s\n",descriptor,value);

    if (strcmp (descriptor, "Lattice_Dimensionality")==0){
      sscanf(value,"%d",&Lattice_Dimensionality);
      if ((Lattice_Dimensionality!=1) &&
          (Lattice_Dimensionality !=2) &&
          (Lattice_Dimensionality !=3)){
          printf("ERROR. Out of range Lattice_Dimensionality. Admissible values are 1, 2 3.\n");
          exit(0);
        }
    }

    else if (strcmp (descriptor, "L")==0){
      sscanf(value,"%d",&L[n_box_sides]);
      n_box_sides++;
    }
  }

  /* Check that the number of specified box sides matches the declared dimensionality */

   if (Lattice_Dimensionality!= n_box_sides){
    printf("Error. Number of specified lattice sides (%d) does not match declared lattice dimensionality (%d)\n",n_box_sides,Lattice_Dimensionality);
    exit(1);
   }


  printf("Finished reading lattice parameters.\n\n");
  printf("Summary of lattice parameters:\n");
  printf("\tLattice_Dimensionality: %d\n",Lattice_Dimensionality);
  printf("\tL: "); for (i=0; i < Lattice_Dimensionality; i++) printf(" %d",L[i]); printf("\n");

  printf("\n\n");

}

/****************************************/

void read_input_temperatures(char *filename, struct System_properties system[]){

  int i;
  FILE *fp;

  printf("Reading %d Replica temperatures from input file %s\n", N_replicas, filename);
  fp=fopen(filename,"r");
  for(i=0; i < N_replicas; i++){
    fscanf(fp,"%lf",&(system[i].temperature));
  }
  fclose(fp);

  for(i=0; i < N_replicas; i++){
    printf("Temperature[%d]: %lf\n",i,system[i].temperature);
  }
  printf("\n\n");

}

/****************************************/


void  dump_observables(struct Lattice *lattice,  struct System_properties *system, int Output_trj, int Output_energy, int Output_spin_config){

int i, j, k, site0, site1;
double coord[3], scale;
int discretization;
FILE *fp;

if (Output_trj==1){
  discretization=10;
  scale=10.0;
  fp=fopen(system->trj_filename,"a");
  fprintf(fp,"%d\n\n",discretization*lattice->n_sites);
  for(i=0; i < system->n_spins; i++) {
  //  printf("spin[%d]= %d\n",i,system->spin_array[i]);
    if (system->spin_array[i]==1){
      site0= lattice->bond[i][0];
      site1= lattice->bond[i][1];
      for(j=0; j < discretization; j++){
        for(k=0; k < 3; k++){
          coord[k]=lattice->site[site0][k] + j*1.0*(lattice->site[site1][k]- lattice->site[site0][k])/discretization;
        }
        fprintf(fp,"C %lf %lf %lf\n",scale*coord[0],scale*coord[1],scale*coord[2]);
      }
    }
  }
  fclose(fp);
}


if (Output_energy==1){
  fp=fopen(system->energy_filename,"a");
  fprintf(fp,"%d\n",system->energy);
  fclose(fp);
}

if (Output_spin_config==1){
  fp=fopen(system->configs_filename,"a");
  fprintf(fp,"%d    ",system->energy);
  for(i=0; i < system->n_spins; i++) {
    fprintf(fp,"%d ",system->spin_array[i]);
  }
  fprintf(fp,"\n");
  fclose(fp);
}

}


