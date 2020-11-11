// show_conf.c
// print out the location of obstacles and particles at given config #

#include "setup.h"

char directory[100];

int main(int argc, char **argv) {  

  int iblo = -1;
  switch (argc){
    case 3:
      sscanf(argv[2],"%d",&iblo);
      sscanf(argv[1],"%s",directory);
      break;
    default:
      print_and_exit("Usage: %s directory iblo\n", argv[0]);
  }
  
  read_conf(directory, iblo);
  int i;

  for(i=0; i<N; i++)
    printf("%.8g %.8g %d\n",u[i].r[0], u[i].r[1], 0);

  for(i=0; i<M; i++) 
    printf("%.8g %.8g %d\n",obst[i].r[0], obst[i].r[1], 1);
  
  printf("#\tL\t%lf\n",data.L);
  printf("#\tN\t%u\n",data.N);
  printf("#\tnblo\t%u \n",data.nblo);
  printf("#\tmedblo\t%u \n",data.medblo);
  printf("#\tmesfr\t%u \n",data.mesfr);
  printf("#\tzseed\t%llu \n",data.seed);
  printf("#\todiff\t%f \n",odiff);
  printf("#\tR\t%f \n",data.R);
  printf("#\tM\t%u \n",M);
  printf("#\tDim \t%d\n",DIM);

  return 0; //end
}

void free_data() {}
