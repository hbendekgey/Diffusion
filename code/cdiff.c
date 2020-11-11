// cdiff.c
// Simulate the diffusion of a particle in a crowded space
#include "setup.h"

// parameters for random number generation
#define FNORM (5.4210108624275218e-20)// max double such that RAND_MAX*FNORM<1
#define FRANDOM (randcong64()*FNORM) //A random number 0<= r <1
unsigned long long randcong64(void); //The PRNG
int extra_gaussian; // boolean for if there is an extra gaussian value
double gaussian; // leftover N(0,1) value from Box Muller 

int nCells; // number of cells in each row/col of the grid
int totalCells; // total number of cells in the grid, i.e. nCells^Dim
int maxInCell;  // Max number of obstacles allowed in a single cell (crashes if exceeded)
int nn;     // number of neighbors of each cell, counting yourself
int* occupancy; //2d array counting the number of obstacles in each cell
int* neighbors; //3d array (nCells x nCells x nn)containing the neighbors to each cell
Obstacle** grid; //3d array (nCells x nCells x maxInCell) containing pointers to obstacles
Particle* usave; //2d array (N x 1000) containing configurations saved since last block 
int occ_alloced = 0; 
int neighbors_alloced = 0; // booleans to check for mallocs, to free at the end
int grid_alloced = 0;
int usave_alloced = 0;
int nblo;   // Number of simulation blocks
int nconf; // number of configs saved this block
// array of boolean values saying if the particles are slaved by obstacles.
// Because the obstacles only move every OFREQ steps, this value is constant
// across that many iterations, so it is more efficient to store
int covered[MAXN];
// for each particle, contains the squared distance to the nearest obstacle center
// so that if a proposed move is shorter than this, 
// there's no need to check for obstacle overlap
double mindist[MAXN];

void init(double);
void read_data(void);
void init_grid(void);
int getCell(double[]);
void addToCell(Obstacle *, int);
int getGridIndex(int, int);
double getNearestObst(double[]);
double getGaussian();
void moveParticle(int);
void moveObstacle(int);
void measureDR2(double *, int);
void save_conf(double *, int, double);

int main(int argc, char **argv) {  
  int iblo = 0, imes, it, i; // by default, start from beginning
  double d, odensity; // how many obstacle centers per unit area 
  nblo = 0; // default to running 0 blocks

  switch (argc){
    case 4:
      sscanf(argv[3],"%d",&nblo);
    case 3:
      sscanf(argv[1],"%lf",&odiff);
      sscanf(argv[2],"%lf",&odensity);
      break;
    default:
      print_and_exit("Usage: %s odiff odensity [nblo]\n", argv[0]);
  }
  //nblo is optional, if we do not specify on the command line, read from input

  // if you see a file named "conf" that contains a block num, 
  // import the configuration and start from there.
  FILE *Fconfig;
  if ((NULL!=(Fconfig=fopen(cfile,"r"))) && (fscanf(Fconfig, "%d", &iblo) == 1)) {
    read_conf(".", iblo);
    init_grid();
    printf("Starting from configuration %d\n", iblo);
    for (i=0; i<N; i++) {
      d = getNearestObst(u[i].r);
      covered[i] = (d <= R2); 
      mindist[i] = d;
    }
  } else {
    printf("Starting from initialization\n");
    init(odensity); //Read input file, initialise PRNG and configuration, etc.
    save_conf(NULL, 0, odensity); // save the initial configuration
  }

  if (NULL != Fconfig) {
    fclose(Fconfig);
  }

  if(nblo)
    data.nblo = nblo;

  //Print parameters to screen
  printf("L\t%lf\n",data.L);
  printf("N\t%u\n",data.N);
  printf("M\t%u\n",M);
  printf("nblo\t%u \n",data.nblo);
  printf("medblo\t%u \n",data.medblo);
  printf("mesfr\t%u \n",data.mesfr);
  printf("seed\t%llu \n",data.seed); 
  printf("odiff\t%f \n",odiff);
  printf("Dim \t%d\n",DIM);

  // import the file of log times to store configs
  FILE *Ftimes;
  if (NULL==(Ftimes=fopen(tfile,"r"))) {
    print_and_exit("Could not open %s", tfile);
  }
  for (i = 0; i < MAXT; i++) {
    fscanf(Ftimes, "%llu", list_times+i);
  }
  fclose(Ftimes);

  // create array of log time configs
  usave = malloc(1000 * N * sizeof(Particle));
  if (usave == NULL) {
    print_and_exit("Unable to allocate memory for usave matrix\n\a");
  } else {
    usave_alloced = 1;
  }

  // array of dr2 measurements
  double measurements[data.medblo];
  //Start simulation
  for(; iblo< data.nblo; iblo++) {
    // to make sure results are perfectly reproducible
    // start fresh w/gaussian generation each block
    extra_gaussian = 0; 
    nconf = 0;
    //Do medblo measurements
    for(imes=0; imes<data.medblo; imes++) {
      for(it=0; it<data.mesfr; it++) { 
        //Do data.mesfr updates without measuring anything
        if(t == list_times[nt]) {
          memcpy(usave+(nconf++)*N, u, (size_t) (sizeof(u[0])*N));
          nt++;
        }
        for(i=0; i<N; i++) {
          if (!covered[i]) {
            moveParticle(i);
          }
        }
        // move obstacles
        if (it % OFREQ == OFREQ - 1) {
          for(i = 0; i < totalCells; i++) {
            occupancy[i] = 0;
          }
          for(i=0; i<M; i++) {
            moveObstacle(i);
          }
          for(i=0; i<N; i++) {
            d = getNearestObst(u[i].r);
            covered[i] = (d <= R2); 
            mindist[i] = d;
          }
        }
        t++;
      }
      if (iblo >= BURNIN) { // don't record dr2 for the first BURNIN blocks
        measureDR2(measurements, imes);
      } else {
        measurements[imes] = 0.;
      }
    }
    save_conf(measurements, iblo+1, odensity);
    //Print some information about the simulation progress
    if (iblo >= BURNIN) {
      printf("Finished block %04d / %04d, dr2 = %1.3lf, log10(dr2/t) = %1.3lf . . .\n", 
            iblo+1, 
            data.nblo, 
            measurements[data.medblo-1], 
            log10(measurements[data.medblo-1]/((iblo - BURNIN + 1) * data.medblo * data.mesfr))); 
    } else {
      printf("Finished block %04d / %04d. Burn-in will take %d blocks\n", iblo+1, data.nblo, BURNIN);
    }
   }

  free_data();
  return 0; 
}

// if we're starting from the beginning:
void init(double odensity) {
  int i, j;
  double d;

  //Read simulation parameters
  read_data();

  //If there is no seed in data, we read a random one
  if(!data.seed)
  {
    FILE *frandom;
    frandom=fopen("/dev/urandom","r");
    fread(&data.seed,(size_t) 8,(size_t) 1,frandom);
    fclose(frandom);
  }
  zseed = data.seed;

  // calculate number of obstacles
  M = floor(odensity * pow(L, DIM));

  for(i = 0; i<M; i++) {
    for (j = 0; j<DIM; j++) {
      obst[i].r[j] = FRANDOM*L;
    }
  }

  init_grid();

  // proposed starting location for the particles
  // we want all particles initialized in the void,
  // so if prop is covered, find a new starting location
  double prop[DIM];
  //Initialise configuration
  for(i=0; i<N; i++) {
    for(j = 0; j < DIM; j++) {
      prop[j] = FRANDOM*L;
    }
    d = getNearestObst(prop);
    if (d <= R2) {
      i--;
      continue;
    }
    covered[i] = 0;
    mindist[i] = d;
    for (j = 0; j < DIM; j++) {
      u[i].r[j] = prop[j];
      u[i].w[j] = 0;
    }
  }
  
  // time parameters for conf_dr2 logging
  t = 0;
  nt = 0;
}

// read the input file
void read_data(void) {

  char coso[1024];
  int j;
  unsigned int * ptdata_int;
  unsigned long long * ptdata_llu;
  double * ptdata_real;
  FILE *Finput;
  Finput=fopen(ifile,"r");
  if (Finput==NULL)
    print_and_exit("I could not open %s\n\a", ifile);

  for (j=0,ptdata_real=&data.L;j<NDATA_DOUBLE;j++)
    {
        fgets(coso,1024,Finput);
        sscanf(coso,"%lf",ptdata_real++);
    }
  for (j=0,ptdata_int=&data.N;j<NDAT_INT;j++)
    {
        fgets(coso,1024,Finput);
        sscanf(coso,"%u",ptdata_int++);
    }
  for (j=0,ptdata_llu=&data.seed;j<NDAT_LLU;j++)
    {
        fgets(coso,1024,Finput);
        sscanf(coso,"%llu",ptdata_llu++);
    }

  fclose(Finput);

  //Shorthands for common variables
  N = data.N;
  L = data.L; 
  R = data.R;
  R2 = pow(R, 2.0);
}

// initialize grid
void init_grid() {
  int i, cell;
  maxInCell = M;
  nCells = ceil(L/R);
  totalCells = (int) pow(nCells, DIM);
  occupancy = malloc(totalCells*sizeof(int));
  if (occupancy == NULL) {
    print_and_exit("Unable to allocate memory for occupancy matrix\n\a");
  } else {
    occ_alloced = 1;
  }
  grid = malloc(totalCells*maxInCell*sizeof(Obstacle*));
  if (grid == NULL) {
    print_and_exit("Unable to allocate memory for grid matrix\n\a");
  } else {
    grid_alloced = 1;
  }
  for(i = 0; i < totalCells; i++) {
    occupancy[i] = 0;
  }
  for (i =0; i < M; i++) {
    addToCell(obst+i,getCell(obst[i].r));
  }
  nn = (int) pow(3.0, DIM) - 1;
  neighbors = malloc(totalCells*nn*sizeof(int));
  if (neighbors == NULL) {
    print_and_exit("Unable to allocate memory for neighbors matrix\n\a");
  } else {
    neighbors_alloced = 1;
  }
  // initialize the array of neighbors
  i = 0; 
  for (cell = 0; cell < totalCells; cell++) {
    #if DIM == 1
      neighbors[i++] = cell==0 ? totalCells-1: cell-1; // left
      neighbors[i++] = cell==totalCells-1 ? 0: cell+1; // right
    #elif DIM == 2
      int j, tempCell;
      for (j = -1; j <= 1; j++) {
        tempCell = cell + j * nCells;
        tempCell = tempCell < 0 ? tempCell + totalCells : tempCell % totalCells;
        neighbors[i++] = tempCell%nCells == 0 ? tempCell+nCells-1 : tempCell-1;
        if (j!= 0)
          neighbors[i++] = tempCell;
        neighbors[i++] = tempCell%nCells == nCells-1 ? tempCell-nCells+1 : tempCell+1;
      }

    #elif DIM == 3
      int j, k, jCell, tempCell, nc2 = nCells*nCells;
      for (j = -1; j <= 1; j++) {
        jCell = cell + j * nc2;
        jCell = jCell < 0 ? jCell + totalCells : jCell % totalCells;
        for (k = -1; k <= 1; k++) {
          if ((jCell % nc2) + k * nCells < 0) {
            tempCell = jCell + nc2 - nCells;
          } else if ((jCell % nc2) + k * nCells >= nc2) {
            tempCell = jCell - nc2 + nCells;
          } else {
            tempCell = jCell + k * nCells;
          }
          neighbors[i++] = tempCell%nCells == 0 ? tempCell+nCells-1 : tempCell-1;
          if (j!= 0 || k != 0)
            neighbors[i++] = tempCell;
          neighbors[i++] = tempCell%nCells == nCells-1 ? tempCell-nCells+1 : tempCell+1;
        }
      }
    #endif
  }
}

// get cell containing o
int getCell(double r[DIM]) {
  int i, cell = 0;
  for (i =0; i < DIM; i++) {
    cell *= nCells;
    cell += floor(r[i]/R);
  } // e.g. o.r[0] * nCells^2 + o.r[1] * nCells + o.r[2] for 3D
  if (cell > totalCells - 1 || cell < 0) {
    print_and_exit("Rounding error with getCell. Improve rounding\n\a");
  }
  return cell;
}

// add obstacle o to the cell
void addToCell(Obstacle *o, int cell) {
  if (occupancy[cell] == maxInCell) {
    print_and_exit("Too much clustering of obstacles; try increasing maxInCell\n\a");
  }
  grid[getGridIndex(cell, occupancy[cell]++)] = o; // increment occupancy and add obstacle
  return;
}

int getGridIndex(int cell, int occ) {
  return cell * maxInCell + occ;
}

// get minimum squared distance to an obstacle (in your cell or neighboring cell).
// note that there might be a closer obstacle that's 2 cells over!!
double getNearestObst(double r[DIM]) {
  int i, j, pcell = getCell(r);
  int *neighbor;
  // default min(d) to the diagonal of a cell
  double dx, d, mind = pow(R,2.) * (DIM+1); 
  // start by checking current cell
  for (i = 0; i < occupancy[pcell]; i++) {
    d = 0.0;
    for (j = 0; j < DIM; j++) {
      dx = fabs(grid[getGridIndex(pcell,i)]->r[j] - r[j]);
      d += pow(dx, 2.0);
    }
    if (d < mind) {
      mind = d;
    }
  }
  // now check all the neighbors
  for (neighbor = neighbors+pcell*nn; neighbor<neighbors+(pcell+1)*nn; neighbor++) {
    for (i = 0; i < occupancy[*neighbor]; i++) {
      d = 0.0;
      for (j = 0; j < DIM; j++) {
        dx = fabs(grid[getGridIndex(*neighbor,i)]->r[j] - r[j]);
        dx = dx < L-dx ? dx : L-dx;
        d += pow(dx, 2.0);
      }
      if (d < mind) {
        mind = d;
      }
    }
  }
  return mind;
}

// Box Muller generates two N(0,1) values at a time. 
// Return one and store the other, unless there's a value already stored
double getGaussian() {
  if (extra_gaussian) {
    extra_gaussian = 0;
    return(gaussian);
  }
  extra_gaussian = 1;
  double r = sqrt(-2 * log(FRANDOM));
  double theta = 2 * M_PI * FRANDOM;
  gaussian = r * cos(theta);
  return(r * sin(theta));
}

// if our current or new position is in a disk, don't move
void moveParticle(int i) {
  double d;
  double prop[DIM];
  // sample a gaussian distance
  double r = PDIFF * getGaussian();
  // then sample a direction to move in
  #if DIM > 2
    double z = 2 * (FRANDOM - 0.5);
    prop[2] = u[i].r[2] + z * r;
  #else 
    double z = 0.0;
  #endif
  #if DIM > 1
    double theta = FRANDOM * M_PI;
    prop[1] = u[i].r[1] + sqrt(1-pow(z,2.0)) * sin(theta) * r;
  #else
    double theta = 0.0;
  #endif 
  prop[0] = u[i].r[0] + sqrt(1-pow(z,2.0)) * cos(theta) * r;  
  int w[DIM];
  int j;
  for (j = 0; j < DIM; j++) {
    w[j] = 0;
    if (prop[j] < 0) {
      prop[j] += L;
      w[j] = -1;
    } else if (prop[j] >= L) {
      prop[j] -= L;
      w[j] = 1;
    }
  }
  // if you're moving less than the min dist to nearest obstacle and are staying
  // in the same cell, don't even check for obstacles
  if (pow(fabs(r)+R, 2) >= mindist[i] || getCell(prop) != getCell(u[i].r)) {
    d = getNearestObst(prop);
    if (d <= R2) {
      return;
    } else {
      mindist[i] = d;
    }
  } else {
    mindist[i] = 0.;
  }
  for (j = 0; j < DIM; j++) {
    u[i].r[j] = prop[j];
    u[i].w[j] += w[j];
  }
}

void moveObstacle(int i) {
  // sample a gaussian distance
  double r = odiff * getGaussian();
  // then sample a direction to move in
  #if DIM > 2
    double z = 2 * (FRANDOM - 0.5);
    obst[i].r[2] += z * r;
  #else 
    double z = 0.0;
  #endif
  #if DIM > 1
    double theta = FRANDOM * M_PI;
    obst[i].r[1] += sqrt(1-pow(z,2.0)) * sin(theta) * r;
  #else
    double theta = 0.0;
  #endif 
  obst[i].r[0] += sqrt(1-pow(z,2.0)) * cos(theta) * r;

  int j;
  for (j = 0; j < DIM; j++) {
    if (obst[i].r[j] < 0) {
      obst[i].r[j] += L;
    } else if (obst[i].r[j] >= L) {
      obst[i].r[j] -= L;
    }
  }
  addToCell(obst + i, getCell(obst[i].r));
}

void measureDR2(double *measurements, int iblo) {
  double dr2 = 0.0;
  int i;
    for(i=0; i<N; i++)
    {
      dr2 += computeSD(u0[i],u[i])/N;
    }
  measurements[iblo] = dr2;
}

void save_conf(double *measurements, int jblo, double odensity) {

  char hostname[100];
  pid_t pid;
  //We want to save a complete snapshot
  //This includes the particles' position,
  //but also the state of the PRNG
  //and also simulation parameters for
  //consistency checks

  FILE *Fconfig;
  char name[1024];
  gethostname(hostname,99);
  pid=getpid();
  sprintf(name,"conf_%06d",jblo);
  if(NULL==(Fconfig=fopen(name,"wb")))
    print_and_exit("I could not open %s\n\a",name);

  //We save the basic parameters
  fwrite(hostname,sizeof(char),99,Fconfig);
  fwrite(&pid,sizeof(pid_t),1,Fconfig);
  fwrite(&data.L,sizeof(data),1,Fconfig);
  // command line params
  fwrite(&odiff, sizeof(double),1,Fconfig);
  fwrite(&M, sizeof(int),1,Fconfig);
  //Current state of the PRNG
  fwrite(&zseed, sizeof(zseed), 1, Fconfig);
  fwrite(&t, sizeof(t), 1, Fconfig);
  fwrite(&nt, sizeof(nt), 1, Fconfig);
  //Particle positions
  fwrite(u, sizeof(Particle), (size_t) N, Fconfig);
  fwrite(obst, sizeof(Obstacle), (size_t) M, Fconfig);

  fclose(Fconfig);

  int i, j;
  FILE *Fmeasure;
  FILE *Ftimes;
  if (jblo == 0) {
    // opening the file with wb will overwrite whatever is currently there,
    // so that later when we open the file to append to it, we're not adding to old stuff
    if(NULL==(Fmeasure=fopen(mfile,"wb")))
      print_and_exit("I could not open %s\n\a", mfile);
    fclose(Fmeasure);
    if (NULL==(Ftimes=fopen(cdfile,"wb")))
      print_and_exit("Could not open %s", cdfile);
    fclose(Ftimes);
  } else {
    // for non-initial save_confs, write the log configs that were stored
    if (NULL==(Ftimes=fopen(cdfile,"ab"))) {
      print_and_exit("Could not open %s", cdfile);
    }
    for (i = 0; i < nconf; i++) {
      fwrite(list_times+nt-nconf+i, sizeof(unsigned long long), 1, Ftimes);
      fwrite(usave+i*N, sizeof(Particle), N, Ftimes);
    }
    fclose(Ftimes);
  }

  if (jblo > BURNIN) {
    // for all post-burnin save_configs, write measurements to dr2.dat
    // and write log-time configs to conf_dr2
    if(NULL==(Fmeasure=fopen(mfile,"a")))
      print_and_exit("I could not open %s\n\a", mfile);
    // each time we write to mfile, the first line's ddr2 must be calculated differently
    double ddr2 = (measurements[0] - last_measure)/data.mesfr;
    fprintf(Fmeasure, "%1.3lf %1.3lf %lf %lf %lf %d\n", measurements[0], ddr2, odensity, odiff,
                                               L, (jblo-BURNIN-1)*data.medblo*data.mesfr+data.mesfr);
    for (i = 1; i < data.medblo; i++) {
      ddr2 = (measurements[i] - measurements[i-1])/data.mesfr;
      fprintf(Fmeasure, "%1.3lf %1.3lf %lf %lf %lf %d\n", measurements[i], ddr2, odensity, odiff,
                                               L, (jblo-BURNIN-1)*data.medblo*data.mesfr+(i+1)*data.mesfr);
    }
    last_measure = measurements[data.medblo - 1];
    fclose(Fmeasure);
  } else {
    // during the burnin period, or for initial configs, update/write u0
    last_measure = 0;

    // update u0
    for (i=0; i<N; i++) {
      for (j=0; j<DIM; j++) {
        u0[i].r[j] = u[i].r[j];
        u0[i].w[j] = u[i].w[j];
      }
    }

    // also write the starting configuration to a file, for calculation of dr2
    FILE * Finit;
    if(NULL==(Finit=fopen(cifile,"wb")))
      print_and_exit("I could not open %s\n\a",cifile);
    fwrite(u0, sizeof(Particle), (size_t) N, Finit);
    fclose(Finit);
  }

  // write file with hardlink to most recent conf
  if(NULL==(Fconfig=fopen(cfile,"wb")))
    print_and_exit("I could not open %s\n\a", cfile);
  fprintf(Fconfig, "%d\n", jblo);
  fclose(Fconfig);
}

//A simple PRNG (we can change it later)
unsigned long long randcong64(void) {
  return zseed=zseed*3202034522624059733LLU+1;
}

void free_data() {
  if(occ_alloced) {
    free(occupancy);
  }
  if(neighbors_alloced) {
    free(neighbors);
  }
  if(grid_alloced) {
    free(grid);
  }
  if(usave_alloced) {
    free(usave);
  }
}
