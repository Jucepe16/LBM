//Incluir algunas librerias que se necesitan
#include <math.h>
#include <stdio.h>
//#include "seconds.h"
#include <complex.h>
#include <iostream>

#include <fstream>  // file streams
#include <sstream>  // string streams
#include <cstdlib>  // standard library

using namespace std;

int maximo(int,int,int,int,int,int);

int main(){

int a=2;
int b=5;
int c=74;
int d=-3;
int e=22;

int max=maximo(5, a, b, c , d, e);
printf("maximo de los numeros %d, %d, %d, %d, %d: %d\n", a,b,c,d,e,max);
//cout<<"maximo:"<<max<<endl;

return 0;
}

int maximo(int n, int r1, int r2, int r3, int r4, int r5){
  int max=0;
  int a[n]={r1, r2, r3, r4,  r5};

  for (int i=0; i<n;i++){
    max=((max> a[i])? max: a[i]);
  }
  return max;
}
