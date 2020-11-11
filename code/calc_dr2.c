// calc_dr2.c
// reads conf_dr2, which has system configurations stored at time steps
// enumerated in powt.txt, and calculates the mean squared displacement
// by averaging across a sliding time window;
// e.g. dr2 for dt=1 is calculated by comparing t0 w/ t1, t1000000 w/ t1000001,
// t2000000 w/ t2000001, etc.
#include "setup.h"

#define MAXTINDEX 3000
#define MAXPOWT 27

double displacements[MAXPOWT]; // log times imported from tfile
int ndisplacements[MAXPOWT]; // number of sliding windows averaged to get the above array
Particle confs[1000000];
unsigned long long utimes[MAXTINDEX];

void get_log_displacements();
void calc_dr2(int, int, int);

char directory[100];

int main(int argc, char **argv) {  
  switch (argc) {
    case 2:
      sscanf(argv[1],"%s",directory);
      break;
    default:
      print_and_exit("Usage: %s directory\n", argv[0]);
  }
  
  read_conf(directory, 0);
  get_log_displacements();

  int i;
  for (i = 0; i < MAXPOWT && ndisplacements[i] > 0; i++) {
    printf("%d %.8g %d %lf\n",(int) pow(2.,i), displacements[i]/ndisplacements[i], ndisplacements[i], odiff);
  }

  return 0; //end
}

void get_log_displacements() {

  int i, j, pt;
  for (i = 0; i < MAXPOWT; i++) {
    displacements[i] = 0.0;
    ndisplacements[i] = 0;
  }

  // read all the configs
  FILE *Ftimes;
  char name[1024];
  sprintf(name,"%s/%s",directory, cdfile);
  if (NULL==(Ftimes=fopen(name,"rb"))) {
    print_and_exit("Could not open %s", name);
  }
  for (i = 0; i < MAXTINDEX; i++) {
    if (fread(utimes+i, sizeof(unsigned long long), 1, Ftimes)) {
      fread(confs+i*N, sizeof(Particle), N, Ftimes);
    } else {
    break;
    }
  }

  int maxT = i;
  fclose(Ftimes);
  for (i=0; i < maxT; i++) {
    if (utimes[i] < BURNIN * data.mesfr * data.medblo) {
      continue;
    }
    if (utimes[i] % 1000000 == 0) { // if we are at a "starting time"
      pt = 0;
      for (j = i+1; j < maxT; j++) {
        if (utimes[i] + (unsigned long long) pow(2.,pt) == utimes[j]) {
          // if the time difference is a power of 2
          calc_dr2(i,j,pt);
          pt++;
        }
      }
    }
  }
}

// calculate square displacement between configs i and j, with dt = 2^pt,
// and update the global arrays
void calc_dr2(int i, int j, int pt) {
  int k;
  for (k = 0; k < N; k++) {
    displacements[pt] += computeSD(confs[i*N+k],confs[j*N+k])/N;
  }
  ndisplacements[pt]++;
}

void free_data() {}
