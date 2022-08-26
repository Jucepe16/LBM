//Incluir algunas librerias que se necesitan
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <iostream>

#include <fstream>  // file streams
#include <sstream>  // string streams
#include <cstdlib>  // standard library

using namespace std;

int nx=2;
int ny=2;

double **a, **b, **auxa, **auxb;

double sum_a2=0.;
double sum_b2=0.;
double sum_dif_a2=0.;
double sum_dif_b2=0.;
double error_a=0.;
double error_b=0.;

int main(){

  //inicializar variables
  a=new double *[nx];
  b=new double *[nx];
  auxa=new double *[nx];
  auxb=new double *[nx];

  for(int i=0; i<ny; i++){
    a[i]=new double [ny];
    b[i]=new double [ny];
    auxa[i]=new double [ny];
    auxb[i]=new double [ny];
  }

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      a[i][j]=0.;
      b[i][j]=0.;
      auxa[i][j]=0.;
      auxb[i][j]=0.;
    }
  }

  //ciclo temporal
  for (int t=0; t<4; t++){

    printf("T=%d\n", t);

    //Asignar valores del paso t
    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
        a[i][j]=(2.*j-0.5)*(t+1);
        b[i][j]=(3.*i-2*j+2.)*(t+1);

        printf("a[%d][%d]= %lf\n",i,j,a[i][j]);
        printf("b[%d][%d]= %lf\n",i,j,b[i][j]);
        printf("auxa[%d][%d]= %lf\n",i,j,auxa[i][j]);
        printf("auxb[%d][%d]= %lf\n",i,j,auxb[i][j]);

        sum_a2+=a[i][j]*a[i][j];
        sum_b2+=b[i][j]*b[i][j];
        sum_dif_a2+=(a[i][j]-auxa[i][j])*(a[i][j]-auxa[i][j]);
        sum_dif_b2+=(b[i][j]-auxb[i][j])*(b[i][j]-auxb[i][j]);

      }
    }

    error_a=sqrt(sum_dif_a2/sum_a2);
    error_b=sqrt(sum_dif_b2/sum_b2);
    printf("Error en a: %lf\n",error_a);
    printf("Error en b: %lf\n",error_b);

    //Asignamos los valores del paso t a los auxiliares para usarlos en el paso t+1
    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
        auxa[i][j]=a[i][j];
        auxb[i][j]=b[i][j];
      }
    }

    sum_a2=0.;
    sum_b2=0.;
    sum_dif_a2=0.;
    sum_dif_b2=0.;

  }

  return 0;
}
