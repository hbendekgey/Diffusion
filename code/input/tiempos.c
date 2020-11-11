#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>

// necessary to make the powt.txt file for logging diffusivity at logarithmic intervals

int main(int argc, char **argv)
{
  int nt,j;
  unsigned long long t;
  unsigned long long lista_t[150];
  unsigned long long t1;

  nt=0;
  lista_t[nt]=0;
  for(j=0;j<=36;j++){
    t=rint(pow(2.,j));
    if(t>lista_t[nt]){
      nt++;
      lista_t[nt]=t;
    }
  }

  for(j=0;j<10000;j++)
  {
    t1 =j*1e6;
    for(nt=0;nt<36;nt++){
      t=t1+lista_t[nt];
	printf("%llu\n",t);
    }
  }
  return 0;
}
