// setup.h
// Assumptions: mesfr % OFREQ == 0, L % R == 0;
// setup file containing global variables used for both simulation and 
// reading configurations/calculation of displacements

#include <stdio.h> 
#include <stdlib.h>  
#include <math.h>   
#include <unistd.h>
#include <stdarg.h>
#include <string.h>

// if Dimension is not specified, we're in 2D
#ifndef DIM
    #define DIM 2
#elif DIM < 1 || DIM > 3
    #error Dimension must be set to 1, 2, or 3 
#endif

// hard-coded values
#define MAXN 10000 //Maximum number of random walkers
#define MAXM 100000 //Maximum number of obstacles
#define PDIFF 1.0 // particle step size, the primary length scale
#define OFREQ 10 // frequency and which obstacles move
#define MAXT 360000
#define BURNIN 10 // number of blocks to throw out at the start
unsigned long long zseed; //Current seed for the PRNG
void print_and_exit(char *, ...); //Print out an error message and exits

// Data type for particles. r is the location vector and w is the winding number vector
typedef struct{double r[DIM]; int w[DIM];} Particle;
typedef struct{double r[DIM];} Obstacle;

Obstacle obst[MAXM]; //Array of obstacles
Particle u[MAXN]; //Array of particles
Particle u0[MAXN]; //Particle starting configurations
double L;   // Size of the simulation box
double R;   // radius of obstacles.
double R2; // radius squared
int N;      // Number of particles
int M;      // NUmber of obstacles
double last_measure; // the most recent value of dr2, to calculate ddr2  
double odiff; // obstacle step size 
char *mfile = "dr2.dat"; // file to write measurements to
char *cfile = "conf"; // file with hard link to most recent config
char *cifile = "conf_init"; // file with particle starting configuration and neighbors array
char *cdfile = "conf_dr2"; // file with log time configurations for dr2 calculation
char *ifile = "input"; // file with input parameters
char *tfile = "powt.txt"; // file with list of log times
unsigned long long t; // time from the start of the second block; used to store configs at logarithmic intervals
int nt; // current index in the powt array
unsigned long long list_times[MAXT]; // log times imported from tfile

//Data structure with system parameters
typedef struct {
  double L,          // size of the simulation box
    R;               // Radius of obstacles
 unsigned int N,     //Number of particles
    nblo,            // number of data block to be generated per magnetization
    medblo,          // number of saved measurements per block 
    mesfr;           // frecuency of measurements
 unsigned long long
    seed;            // Initial seed
} s_data;
#define NDATA_DOUBLE 2 //Number of double-type data  
#define NDAT_INT 4 //Number of integer-type data in structre of type data
#define NDAT_LLU 1 //Number of unsigned long long in  data

s_data data;

void free_data(void);
double computeSD(Particle, Particle);

// if we're continuing from a previous run:
void read_conf(char *directory, int jblo) {

  char hostname[100];
  pid_t pid;

  FILE *Fconfig;
  char name[1024];
  sprintf(name,"%s/conf_%06d",directory,jblo);
  if(NULL==(Fconfig=fopen(name,"rb")))
    print_and_exit("I could not open %s\n\a",name);

  //We read the basic parameters
  fread(hostname,sizeof(char),99,Fconfig);
  fread(&pid,sizeof(pid_t),1,Fconfig);
  fread(&data.L,sizeof(data),1,Fconfig);
  // command line params
  fread(&odiff, sizeof(double),1,Fconfig);
  fread(&M, sizeof(int),1,Fconfig);
  //Current state of the PRNG
  fread(&zseed, sizeof(zseed), 1, Fconfig);
  fread(&t, sizeof(t), 1, Fconfig);
  fread(&nt, sizeof(nt), 1, Fconfig);
  N = data.N;
  L = data.L;
  R = data.R;
  R2 = pow(R, 2.0);
  //Particle positions
  fread(u, sizeof(Particle), (size_t) N, Fconfig);
  fread(obst, sizeof(Obstacle), (size_t) M, Fconfig);

  fclose(Fconfig);

  // read initial configuration for calculating ddr2
  FILE * Finit;
  sprintf(name,"%s/%s",directory,cifile);
  if(NULL==(Finit=fopen(name,"rb")))
    print_and_exit("I could not open %s\n\a",cifile);
  fread(u0, sizeof(Particle), (size_t) N, Finit);
  fclose(Finit);

  // initialize with last value of dr2 so we can calculate ddr2
  last_measure = 0.0;
  int i;
  for(i=0; i<N; i++) {
    last_measure += computeSD(u0[i],u[i])/N;
  }

}

double computeSD(Particle v0, Particle v) {
  int j;
  double d = 0.0;
  for (j = 0; j < DIM; j++) {
    d += pow(v.r[j] - v0.r[j] + (v.w[j] - v0.w[j]) * L, 2.0);
  }
  return d;
}

void print_and_exit(char *format, ...) {
    free_data();
    va_list list;

    va_start(list,format);
    vprintf(format,list);
    va_end(list);
    exit(1);
}
