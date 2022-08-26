//Simulaciòn de una impureza circular en una red 2D clasica (no relativista)
//aqui ya voy a meter diferentes tipos de condiciones de frontera (variable "BC"):
//1) BB
//2) Equilibriuem Scheme
//3) Non-Equilibrium Extrapolation Method (NEEM)
//4) Non-Equilibrium BB (NEBB)

//Incluir algunas librerias que se necesitan
#include <math.h>
#include <stdio.h>
#include "seconds.h"
#include <complex.h>
#include <iostream>

#include <fstream>  // file streams
#include <sstream>  // string streams
#include <cstdlib>  // standard library

using namespace std;

#define SQ(x) (x * x) // square function; replaces SQ(x) by ((x) * (x)) in the code*/
double PI= atan2(1, 1) * 4;
typedef double _Complex cplx;

//Definiciòn algunos paràmetros de relevancia

int Bc=4;                            //BCs type
int form=1;                         //0= circle; 1= sqare; 2= none.

const int Nsteps=50000;                      // Numero de pasos temporales
const int scale = 50;                     // fija el tamaño de la simulaciòn
const int Npop = 9;                        // Numero de direcciones
const int Nx=50*scale;                      // largo del canal
const int Ny=8*scale;                        // ancho del canal
const double cl=1.0;
const double tau_0=0.1;
const double u_max=.1;                  // velocidad maxima
double c= .5;
double tau= tau_0 + c;              // relaxation time
double omega =1./tau;
double nu= (tau-.5)/3.;                         // viscosidad cinematica de sizalle

// variables para la mascara
const int x_0=Ny/2;
const int y_0=Ny/2;
const int r=25;

const double L=2.*(double)r;
double Re= (double)L*u_max/nu;            // numero de Reynolds

int **type;
int fronteras=0;
double epsilon=1.43;
double epsilon_sq=0.1;
double fD_x=0.;
double fD_y=0.;
double j_x=0.;

//Paràmetros de la Red
int cx[9]={ 0, 1, 0, -1, 0, 1, -1, -1, 1 };
int cy[9]={ 0, 0, 1, 0, -1, 1, 1, -1, -1 };
double w[9]={ 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};

//primero se declara (declararemos todas las funciones y despues se definiran)
void init();
void equilibrium(double,double,double);
void LBM();
void momenta();

//funcion para calcular la Drag force
double drag(double, double);         //(density_c, velocity_x, //nu)

//funciones de FFT
void main_fft(cplx*, int);
void _fft(cplx*, cplx*, int, int);
void fft(cplx*, int);
void show(const char*, cplx*);

//funciones de las MASCARAS
void init_type();
void mask(int, int, int);       //circular mask(x_0,y_0,r)
void mask_sq(int, int, int);    //squared mask(x_0,y_0,r)
void no_mask();

//direcciones inversas
int inv(int);

//declaramos 3 poblaciones que ocupamos: actual(f), propagada(f_prop) y equilibrio(feq)
double ***f, ***f_prop;                 // populations (old/propagated) (*** array 3D)
double **density_c;                     // fluid density (** array 2D)
double **vel_x;                         // fluid velocity (x-component)
double **vel_y;                         // fluid velocity (y-component)
double f_eq[Npop], f_init[Npop];
double F_D[Nsteps/2];                   //drag force
double J_x[Nsteps];                     //current in x-direction
cplx buf[Nsteps/2];                     //calculo de FFT
cplx buf1[Nsteps];

//++++++++++++++++++++++++++++++Metodo main() ++++++++++++++++++++++++++++++++++++++++++++
int main(){

//Archivos con los datos para diversas gràficas
  FILE* perfil= fopen("perfil.dat", "w+");
  FILE* datos = fopen("flow_x.dat", "w+");
  FILE* datos1 = fopen("flow_y.dat", "w+");
  FILE* datos2 = fopen("flow_abs.dat", "w+");
  FILE* FFT_re = fopen("fft_re.dat", "w+");
  FILE* FFT_im = fopen("fft_im.dat", "w+");
  FILE* FFT_abs = fopen("fft_abs.dat", "w+");
  FILE* FFT_current_re = fopen("fft1_re.dat", "w+");
  FILE* FFT_current_im = fopen("fft1_im.dat", "w+");
  FILE* FFT_current_abs = fopen("fft1_abs.dat", "w+");
  FILE* obstacle = fopen("mask.dat","w+");
  FILE* drag = fopen("drag_t.dat", "w+");
  FILE* curr = fopen("curr_t.dat", "w+");

  if(datos==NULL){
    printf("Error en la creaciòn del archivo\n\n");
  }
  else{
    cout<<"INICIA EL PROGRAMA LBM"<<endl;
    printf("Simulating electronic flow in graphene.\n");
    printf("      domain size: %dx%d\n",Nx,Ny);
    printf("               nu: %lf\n",nu);
    printf("              tau: %lf\n",tau);
    printf("            u_max: %lf\n",u_max);
    //printf("                T: %lf\n",T);
    printf("        timesteps: %d\n",Nsteps);
    printf("        Re: %lf\n",(u_max*L/nu));
    //printf("       save every: %u\n",NSAVE);
    //printf("    message every: %u\n",NMSG);
    printf("\n");

    //double bytesPerMiB = 1024.0*1024.0;
    //double bytesPerGiB = 1024.0*1024.0*1024.0;

    init(); //inicializa: rho_c=1, v_x=0.002, v_y=0, f =f_prop=f_eq

    //mask
    if(form==0){
      printf("Se coloca una barrera circular en la posiciòn (%d, %d) de radio %d.\n", x_0, y_0, r);

      mask(x_0, y_0, r);
      for(int X=0; X<Nx; X++){
        for(int Y=0; Y<Ny; Y++){
          //if(type[X][Y]==0 || type[X][Y]==2){
            fprintf(obstacle, "%d, %d, %d\n", X,Y,type[X][Y]);
            //printf("(X,Y)= (%d,%d). SQ(X-x_0)+SQ(Y-y_0)= %d\n", X,Y,(X-x_0)*(X-x_0)+(X-x_0)*(Y-y_0));
            //printf("(X,Y)= (%d,%d). SQ(X-x_0)+SQ(Y-y_0)= %d\n", X,Y,SQ((X-x_0))+SQ((Y-y_0)));
          //}
        }
      }
    }

    else if(form==1){
      printf("Se coloca una barrera cuadrada en la posiciòn (%d, %d) de radio %d.\n", x_0, y_0, r);
      mask_sq(x_0, y_0, r);
      for(int X=0; X<Nx; X++){
        for(int Y=0; Y<Ny; Y++){
          //if(type[X][Y]==0 || type[X][Y]==2){
            fprintf(obstacle, "%d, %d, %d\n", X,Y,type[X][Y]);
            //printf("(X,Y)= (%d,%d). SQ(X-x_0)+SQ(Y-y_0)= %d\n", X,Y,(X-x_0)*(X-x_0)+(X-x_0)*(Y-y_0));
            //printf("(X,Y)= (%d,%d). SQ(X-x_0)+SQ(Y-y_0)= %d\n", X,Y,SQ((X-x_0))+SQ((Y-y_0)));
          //}
        }
      }
    }

    else{
      no_mask();
      for(int X=0; X<Nx;X++){
        for(int Y =0; Y<Ny; Y++){
          fprintf(obstacle, "%d, %d, %d\n", X,Y,type[X][Y]);
        }
      }
      printf("No barrier set.");
    }

    double start = seconds();  //calculo del tiempo de ejecucion

    /*//Imprimir los valores iniciales de la velocidad en el archivo .dat
    for(int X=0; X<Nx;X++){
      for(int Y=0; Y<Ny;Y++){
      //fprintf(datos,"%d %lf\n", Y, vel_x[-1+Nx/2][Y]);
      fprintf(datos,"%d %d %lf\n", X, Y, vel_x[X][Y]);
      }
    }*/

    for(int t=0;t<Nsteps;t++){                        //Nsteps pasos temporales

       LBM();
       //printf("paso: %d, T_0/T: %lf\n",t+1, T_0/T);
       //printf("paso: %d, Tau: %lf\n",t+1, tau);
       //printf("paso: %d, T: %lf, E: %lf\n",t+1, T, density_e[Nx-1][Ny/2]);

       /*//fprintf(datos,"#t= %d\n", t);
       //Imprimir los valores finales de la velocidad en el archivo .dat
       for(int X=0; X<Nx;X++){*/

      // }
       //fprintf(datos,"\n \n");


       //imprimir los valores de la drag force en el arreglo buf[t]
       if((t%2)==0){
       F_D[t]=fD_x;
       //J_x[t]=j_x;
       fprintf(drag, "%d %lf \n", t, F_D[t]);
       fD_x=0.;
     }
     J_x[t]=j_x;
     fprintf(curr, "%d %lf \n", t, J_x[t]);
     j_x=0.;
    }

    //Imprimir los valores finales de la velocidad en el archivo .dat
    //fprintf(datos2,"#X Y Z \n");
    for(int X=0; X<Nx;X++){
      for(int Y=0; Y<Ny;Y++){
      fprintf(perfil,"%d %lf\n", Y, vel_x[-1+Nx/2][Y]);
      fprintf(datos,"%d %d %lf\n", X, Y, vel_x[X][Y]);
      fprintf(datos1,"%d %d %lf\n", X, Y, vel_y[X][Y]);
      fprintf(datos2,"  %d  %d  %lf \n", X, Y, sqrt(SQ(vel_y[X][Y])+SQ(vel_x[X][Y])));
      }
      fprintf(datos,"\n");
      fprintf(datos1,"\n");
      fprintf(datos2,"\n");
    }
    fprintf(datos,"\n");
    fprintf(datos1,"\n");
    fprintf(datos2,"\n");

    for(int t=0; t<Nsteps/2;t++){
      buf[t]=F_D[t];
      //buf1[t]=J_x[t];
    }
    for(int t=0; t<Nsteps;t++){
      //buf[t]=F_D[t];
      buf1[t]=J_x[t];
    }
    //show("Data: ", buf);
    fft(buf, Nsteps/2);
    fft(buf1, Nsteps);
    //show("\nFFT : ", buf);

    for(int t=0; t<Nsteps;t++){
      fprintf(FFT_re,"%lf %lf\n", (double)t/((double)Nsteps-1.), creal(buf[t]));
      fprintf(FFT_im,"%lf %lf\n", (double)t/((double)Nsteps-1.), cimag(buf[t]));
      fprintf(FFT_abs,"%lf %lf\n", (double)t/((double)Nsteps-1.), cabs(buf[t]));
      fprintf(FFT_current_re,"%lf %lf\n", (double)t/((double)Nsteps-1.), creal(buf1[t]));
      fprintf(FFT_current_im,"%lf %lf\n", (double)t/((double)Nsteps-1.), cimag(buf1[t]));
      fprintf(FFT_current_abs,"%lf %lf\n", (double)t/((double)Nsteps-1.), cabs(buf1[t]));
    }

    double end = seconds();
    double runtime= end-start;
    printf("Tiempo de còmputo: %lf segundos.\n", runtime);

    //printf("Numero de puntos frontera: %d\n", fronteras);

    return 0;

    fclose(datos);
    fclose(datos1);
    fclose(datos2);
    fclose(FFT_re);
    fclose(FFT_im);
    fclose(FFT_abs);
    fclose(FFT_current_re);
    fclose(FFT_current_im);
    fclose(FFT_current_abs);
    fclose(curr);
    fclose(drag);

    //de la mask
    fclose(obstacle);

    free(f);  free(f_prop); free(f_eq);
    free(density_c);  free(vel_x); free(vel_y);
    free(f_init);
    free(F_D); free(buf); free(J_x); free(buf1);
  }
}

//***********************************************************************************
//**********init(): ESTA FUNCION ALOJA MEMORIA E INICIALIZA LAS VARIABLES************
//***********************************************************************************
void init(){
    init_type();           //inicializa los tipos de celdas todas iguales a 0
    if(form==0){
      mask(x_0,y_0,r);     //asigna la mascara con tipos de celdas (fluido, frontera, solido)
    }
    else if (form==1){
      mask_sq(x_0,y_0,r);
    }
    else no_mask();

    //inicializar buf[t] y F_D[t], j_x[t]
    for(int t=0; t<Nsteps; t++){
      buf[t]=0.;
      buf1[t]=0.;
      F_D[t]=0.;
      J_x[t]=0.;
    }

    //alojar memoria para la densidad y las velocidades (vel_x y vel_y)
    //density_e = new double*[Nx];
    density_c = new double*[Nx];
    vel_x = new double*[Nx];
    vel_y = new double*[Nx];

    for(int X = 0; X < Nx; X++) {
      density_c[X] = new double[Ny];
      vel_x[X] = new double[Ny];
      vel_y[X] = new double[Ny];
    }

    // Initialize the fluid density and velocity.
    for(int X = 0; X < Nx; X++) {
      //F_D[X]=0.0;
      for(int Y = 0; Y < Ny; Y++) {
        if(type[X][Y]==2){
          density_c[X][Y] = 1.0;
          vel_x[X][Y] = 0.;
          vel_y[X][Y] = 0.;
        }
        else{
          density_c[X][Y] = 0.0;
          vel_x[X][Y] = 0.0;
          vel_y[X][Y] = 0.0;
        }
      }
    }

     // Allocate memory for the populations (array de [Npop][Nx][Ny])
     f = new double**[Npop];
     f_prop = new double**[Npop];

     for(int c_i = 0; c_i < Npop; c_i++) {
       f[c_i] = new double*[Nx];
       f_prop[c_i] = new double*[Nx];

       for(int X = 0; X < Nx; X++) {
         f[c_i][X] = new double[Ny];
         f_prop[c_i][X] = new double[Ny];

         for(int Y = 0; Y < Ny; Y++) {
           f[c_i][X][Y] = 0.;
           f_prop[c_i][X][Y] = 0.;

         }
       }
     }

     //Initialize populations. Use the equilibrium corresponding to the initialized fluid density and velocity.
     for(int X = 0; X < Nx; X++) {
       for(int Y = 0; Y < Ny; Y++) {
         if(type[X][Y]==2){
           equilibrium(density_c[X][Y], vel_x[X][Y], vel_y[X][Y]);

           for(int c_i = 0; c_i < Npop; c_i++) {
             f_prop[c_i][X][Y] = f_eq[c_i] ;
             f[c_i][X][Y] = f_eq[c_i] ;
             f_init[c_i]=f_eq[c_i];         //guardan los valores iniciales de frontera inlet.
           }
         }
       }
     }

     return;
}

// COMPUTE EQUILIBRIUM

void equilibrium(double den_c, double v_x, double v_y) {

  f_eq[0]=w[0]*den_c*(1.-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));
  f_eq[1]=w[1]*den_c*(1.+3*(( v_x    ))/(SQ(cl))+4.5*SQ(( v_x    ))/SQ(SQ(cl))-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));
  f_eq[2]=w[2]*den_c*(1.+3*((     v_y))/(SQ(cl))+4.5*SQ((     v_y))/SQ(SQ(cl))-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));
  f_eq[3]=w[3]*den_c*(1.+3*((-v_x    ))/(SQ(cl))+4.5*SQ((-v_x    ))/SQ(SQ(cl))-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));
  f_eq[4]=w[4]*den_c*(1.+3*((    -v_y))/(SQ(cl))+4.5*SQ((    -v_y))/SQ(SQ(cl))-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));
  f_eq[5]=w[5]*den_c*(1.+3*(( v_x+v_y))/(SQ(cl))+4.5*SQ(( v_x+v_y))/SQ(SQ(cl))-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));
  f_eq[6]=w[6]*den_c*(1.+3*((-v_x+v_y))/(SQ(cl))+4.5*SQ((-v_x+v_y))/SQ(SQ(cl))-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));
  f_eq[7]=w[7]*den_c*(1.+3*((-v_x-v_y))/(SQ(cl))+4.5*SQ((-v_x-v_y))/SQ(SQ(cl))-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));
  f_eq[8]=w[8]*den_c*(1.+3*(( v_x-v_y))/(SQ(cl))+4.5*SQ(( v_x-v_y))/SQ(SQ(cl))-1.5*(SQ(v_x)+SQ(v_y))/(SQ(cl)));

  return;
}

// **********************************
// COMPUTE FLUID DENSITY AND VELOCITY
// **********************************

// This function computes the fluid density and velocity from the populations.
// The velocity correction due to body force is *not* included here.

void momenta() {
  for(int X = 0; X < Nx; X++) {
    for(int Y = 0; Y < Ny; Y++) {
      if(type[X][Y]==2){
        density_c[X][Y] = (f_prop[0][X][Y] + f_prop[1][X][Y] + f_prop[2][X][Y] + f_prop[3][X][Y] + f_prop[4][X][Y] + f_prop[5][X][Y] + f_prop[6][X][Y] + f_prop[7][X][Y] + f_prop[8][X][Y]);
        vel_x[X][Y] = ((f_prop[1][X][Y] + f_prop[5][X][Y] + f_prop[8][X][Y]) - (f_prop[3][X][Y] + f_prop[6][X][Y] + f_prop[7][X][Y])) / density_c[X][Y];
        vel_y[X][Y] = ((f_prop[2][X][Y] + f_prop[5][X][Y] + f_prop[6][X][Y]) - (f_prop[4][X][Y] + f_prop[7][X][Y] + f_prop[8][X][Y])) / density_c[X][Y];

        if(X==Nx/2 && Y==Ny/2){
          j_x= density_c[X][Y]*vel_x[X][Y];
        }
      }
    }
    //F_D[X]=drag(density_c[X][Ny-2], vel_x[X][Ny-2]);
  }

  return;
}

//***********************LBM()**************************
//Se calculan por separado:colisiòn, inlet boundary conditions, propagaciòn y BCs

void LBM(){

    int newx=0;
    int newy=0;
    fD_x=0.;
    fD_y=0.;

      for(int X = 0; X < Nx; X++) {
        for(int Y = 0; Y < Ny; Y++) {
          if(type[X][Y]==2){

          // Compute equilibrium populations.
          equilibrium(density_c[X][Y], vel_x[X][Y], vel_y[X][Y]);

          // collision step
          f[0][X][Y] = f_prop[0][X][Y] * (1. - omega) + f_eq[0] * omega;
          f[1][X][Y] = f_prop[1][X][Y] * (1. - omega) + f_eq[1] * omega;
          f[2][X][Y] = f_prop[2][X][Y] * (1. - omega) + f_eq[2] * omega;
          f[3][X][Y] = f_prop[3][X][Y] * (1. - omega) + f_eq[3] * omega;
          f[4][X][Y] = f_prop[4][X][Y] * (1. - omega) + f_eq[4] * omega;
          f[5][X][Y] = f_prop[5][X][Y] * (1. - omega) + f_eq[5] * omega;
          f[6][X][Y] = f_prop[6][X][Y] * (1. - omega) + f_eq[6] * omega;
          f[7][X][Y] = f_prop[7][X][Y] * (1. - omega) + f_eq[7] * omega;
          f[8][X][Y] = f_prop[8][X][Y] * (1. - omega) + f_eq[8] * omega;

        }
      }
    }

    //Streaming step
    for(int k=0; k<Npop; k++){
      for(int X = 0; X < Nx; X++) {
        for(int Y = 0; Y < Ny; Y++) {
            if(type[X][Y]==2){
              // Streaming step (Periodic in the Y direction)
              newx=(X+cx[k]+Nx)%Nx; //aqui por que es periodico en x?
              newy=(Y+cy[k]+Ny)%Ny;
              f_prop[k][newx][newy]=f[k][X][Y];
            }
        }
      }
    }

    // Inlet/Outlet BC: PBBC (w/ i=0 and i=NX-1 boundary layers)
    for(int Y=1;Y<Ny-1;Y++){

      //inlet
      //density_c=(f_prop[0][0][Y]+f_prop[1][0][Y]+f_prop[2][0][Y]+f_prop[3][0][Y]+f_prop[4][0][Y]+f_prop[5][0][Y]+f_prop[6][0][Y]+f_prop[7][0][Y]+f_prop[8][0][Y])
      vel_x[0][Y]=0.1;
      vel_y[0][Y]=0.0;
      density_c[0][Y]= (1./(1.-vel_x[0][Y]))*(f_prop[0][0][Y]+f_prop[2][0][Y]+f_prop[4][0][Y]+2.*(f_prop[3][0][Y]+f_prop[6][0][Y]+f_prop[7][0][Y]));
      f_prop[1][0][Y]=f_prop[3][0][Y]+(2./3.)*vel_x[0][Y]* density_c[0][Y];
      f_prop[5][0][Y]=f_prop[7][0][Y]+(1./6.)*vel_x[0][Y]* density_c[0][Y]-0.5*(f_prop[2][0][Y]-f_prop[4][0][Y])+0.5*vel_y[0][Y]* density_c[0][Y];
      f_prop[8][0][Y]=f_prop[6][0][Y]+(1./6.)*vel_x[0][Y]* density_c[0][Y]+0.5*(f_prop[2][0][Y]-f_prop[4][0][Y])-0.5*vel_y[0][Y]* density_c[0][Y];

      //outlet
      //density_c=(f_prop[0][0][Y]+f_prop[1][0][Y]+f_prop[2][0][Y]+f_prop[3][0][Y]+f_prop[4][0][Y]+f_prop[5][0][Y]+f_prop[6][0][Y]+f_prop[7][0][Y]+f_prop[8][0][Y])
      vel_x[Nx-1][Y]=0.1;
      vel_y[Nx-1][Y]=0.0;
      density_c[Nx-1][Y]= (1./(1.+vel_x[Nx-1][Y]))*(f_prop[0][Nx-1][Y]+f_prop[2][Nx-1][Y]+f_prop[4][Nx-1][Y]+2.*(f_prop[1][Nx-1][Y]+f_prop[5][Nx-1][Y]+f_prop[8][Nx-1][Y]));
      f_prop[3][Nx-1][Y]=f_prop[1][Nx-1][Y]-(2./3.)*vel_x[Nx-1][Y]* density_c[Nx-1][Y];
      f_prop[6][Nx-1][Y]=f_prop[8][Nx-1][Y]-(1./6.)*vel_x[Nx-1][Y]* density_c[Nx-1][Y]-0.5*(f_prop[2][Nx-1][Y]-f_prop[4][Nx-1][Y])+0.5*vel_y[Nx-1][Y]* density_c[Nx-1][Y];
      f_prop[7][Nx-1][Y]=f_prop[5][Nx-1][Y]-(1./6.)*vel_x[Nx-1][Y]* density_c[Nx-1][Y]+0.5*(f_prop[2][Nx-1][Y]-f_prop[4][Nx-1][Y])-0.5*vel_y[Nx-1][Y]* density_c[Nx-1][Y];

    /*  //inlet in x=0: Missing populations equal to initial conditions
      f_prop[0][0][Y] = f_init[0];
      f_prop[1][0][Y] = f_init[1];
      f_prop[2][0][Y] = f_init[2];
      f_prop[3][0][Y] = f_init[3];
      f_prop[4][0][Y] = f_init[4];
      f_prop[5][0][Y] = f_init[5];
      f_prop[6][0][Y] = f_init[6];
      f_prop[7][0][Y] = f_init[7];
      f_prop[8][0][Y] = f_init[8];

      //equilibrium(density[1][Y], vel_x[1][Y], vel_y[1][Y]);

    //outlet in x=Nx-1: free BCS. Usaremos extrapolaciòn a primer orden
      f_prop[0][Nx-1][Y]= f[0][Nx-2][Y];
      f_prop[1][Nx-1][Y]= f[1][Nx-2][Y];
      f_prop[2][Nx-1][Y]= f[2][Nx-2][Y];
      f_prop[3][Nx-1][Y]= f[3][Nx-2][Y];
      f_prop[4][Nx-1][Y]= f[4][Nx-2][Y];
      f_prop[5][Nx-1][Y]= f[5][Nx-2][Y];
      f_prop[6][Nx-1][Y]= f[6][Nx-2][Y];
      f_prop[7][Nx-1][Y]= f[7][Nx-2][Y];
      f_prop[8][Nx-1][Y]= f[8][Nx-2][Y];
*/
    }

    //corners
    //left-top
    f_prop[1][0][Ny-1]= f[3][0][Ny-1];
    f_prop[4][0][Ny-1]= f[2][0][Ny-1];
    f_prop[8][0][Ny-1]= f[6][0][Ny-1];
    //f_prop[5][0][Ny-1]= f[0][0][Ny-1];
    //f_prop[7][0][Ny-1]= f[0][0][Ny-1];

    //left-bottom
    f_prop[1][0][0]= f[3][0][0];
    f_prop[2][0][0]= f[4][0][0];
    f_prop[5][0][0]= f[7][0][0];
    //f_prop[6][0][0]= f[0][0][0];
    //f_prop[8][0][0]= f[0][0][0];

    //right-top
    f_prop[3][Nx-1][Ny-1]= f[1][Nx-1][Ny-1];
    f_prop[4][Nx-1][Ny-1]= f[2][Nx-1][Ny-1];
    f_prop[7][Nx-1][Ny-1]= f[5][Nx-1][Ny-1];
    //f_prop[6][Nx-1][Ny-1]= f[2][Nx-1][Ny-1];
    //f_prop[8][Nx-1][Ny-1]= f[2][Nx-1][Ny-1];

    //right-bottom
    f_prop[2][Nx-1][0]= f[4][Nx-1][0];
    f_prop[3][Nx-1][0]= f[1][Nx-1][0];
    f_prop[6][Nx-1][0]= f[8][Nx-1][0];
    //f_prop[6][Nx-1][0]= f[6][Nx-1][0];
    //f_prop[7][Nx-1][0]= f[7][Nx-1][0];

  /*    //Streaming step
      for(int k=0; k<Npop; k++){
        for(int X = 0; X < Nx; X++) {
          for(int Y = 0; Y < Ny; Y++) {
              if(type[X][Y]==2){
                // Streaming step (Periodic in the Y direction)
                newx=(X+cx[k]+Nx)%Nx; //aqui por que es periodico en x?
                newy=(Y+cy[k]+Ny)%Ny;
                f_prop[k][newx][newy]=f[k][X][Y];
                }
          }
        }
      }*/

      //Boundary conditions at Top/ybottom
      //1)BB  2)ES  3)NEEM  4)NEBB

      if(Bc==1){
        for(int X=0; X<Nx;X++){
        // Top wall (rest)
          f_prop[4][X][Ny-1]=f[2][X][Ny-1];
          f_prop[7][X][Ny-1]=f[5][X][Ny-1];
          f_prop[8][X][Ny-1]=f[6][X][Ny-1];

          // Bottom wall (rest)
          f_prop[2][X][0]=f[4][X][0];
          f_prop[5][X][0]=f[7][X][0];
          f_prop[6][X][0]=f[8][X][0];

        }
      }

      else{

      //Boundary condition (wet node)

      //Setting macroscopic quantities at boundaries

      //Bottom wall (rest)
      for(int X=0; X<Nx;X++){
        vel_x[X][0]=0.;
        vel_y[X][0]=0.;
        density_c[X][0]=(f_prop[0][X][0]+f_prop[1][X][0]+f_prop[3][X][0]+2.*(f_prop[4][X][0]+f_prop[7][X][0]+f_prop[8][X][0]));

      }

      //Top wall (rest)
      for(int X=0; X<Nx;X++){
        density_c[X][Ny-1]=(f_prop[0][X][Ny-1]+f_prop[1][X][Ny-1]+f_prop[3][X][Ny-1]+2.*(f_prop[2][X][Ny-1]+f_prop[5][X][Ny-1]+f_prop[6][X][Ny-1]));
        vel_x[X][Ny-1]=0.;
        vel_y[X][Ny-1]=0.;
      }

      //Setting populations quantities at boundaries (note: rho=1)

      if(Bc==2){
        for (int k=1; k<Npop; ++k){
          for(int X=0; X<Nx;X++){
            //f_eq[1]   = w[1] * den_c;
            f_prop[k][X][0]   =w[k]* density_c[X][0];
            f_prop[k][X][Ny-1]=w[k]* density_c[X][Ny-1];
          }
        }
      }

      else if(Bc==3){
        for (int k=1; k<Npop; ++k){
          for(int X=0; X<Nx;X++){
            equilibrium(density_c[X][1], vel_x[X][1], vel_y[X][1]);
            //Bottom
            f_prop[k][X][0]   =w[k]* density_c[X][0] +(f_prop[k][X][1]   -f_eq[k]);
            equilibrium(density_c[X][Ny-2], vel_x[X][1], vel_y[X][Ny-2]);
            //Top
            f_prop[k][X][Ny-1]=w[k]* density_c[X][Ny-1]+(f_prop[k][X][Ny-2]-f_eq[k]);
          }
        }
      }

      else{
        for(int X=0; X<Nx;X++){
          //Bottom
          f_prop[2][X][0]   =f_prop[4][X][0]+(2./3.)*vel_y[X][0]*density_c[X][Ny-1];
          f_prop[5][X][0]   =f_prop[7][X][0]+(1./6.)*vel_y[X][0]*density_c[X][Ny-1]-0.5*(f_prop[1][X][0]-f_prop[3][X][0])+0.5*vel_x[X][0]*density_c[X][Ny-1];
          f_prop[6][X][0]   =f_prop[8][X][0]+(1./6.)*vel_y[X][0]*density_c[X][Ny-1]+0.5*(f_prop[1][X][0]-f_prop[3][X][0])-0.5*vel_x[X][0]*density_c[X][Ny-1];

          //Top
          f_prop[4][X][Ny-1]=f_prop[2][X][Ny-1]-(2./3.)*vel_y[X][Ny-1]* density_c[X][Ny-1];
          f_prop[7][X][Ny-1]=f_prop[5][X][Ny-1]-(1./6.)*vel_y[X][Ny-1]* density_c[X][Ny-1]+0.5*(f_prop[1][X][Ny-1]-f_prop[3][X][Ny-1])-0.5*vel_x[X][Ny-1]* density_c[X][Ny-1];
          f_prop[8][X][Ny-1]=f_prop[6][X][Ny-1]-(1./6.)*vel_y[X][Ny-1]* density_c[X][Ny-1]-0.5*(f_prop[1][X][Ny-1]-f_prop[3][X][Ny-1])+0.5*vel_x[X][Ny-1]* density_c[X][Ny-1];

        }
      }
      }

///*
      /*//Streaming step
      for(int k=0; k<Npop; k++){
        for(int X = 0; X < Nx; X++) {
          for(int Y = 0; Y < Ny; Y++) {
              if(type[X][Y]==2){
                // Streaming step (Periodic in the Y direction)
                newx=(X+cx[k]+Nx)%Nx;
                newy=(Y+cy[k]+Ny)%Ny;
                f_prop[k][newx][newy]=f[k][X][Y];
            }
          }
        }
      }
*/

        //BB BCs en el obstàculo sòlido
        for(int k=0; k<Npop; k++){
          for(int X=((x_0-r-2));X<((x_0+r+2));X++){
            for(int Y=((y_0-r-2));Y<((y_0+r+2));Y++){
              if(type[(X+cx[k])][(Y+cy[k])]==1){           //Si es un nodo frontera
                 f_prop[inv(k)][X][Y]=f[k][X][Y];          //hacer la funcion inv(k) que mande la direcciòn
                 fD_x += 2.*f[k][X][Y]*cx[k];              //Drag force
                 fD_y += 2.*f[k][X][Y]*cy[k];              //Lift force
              }
            }
          }
        }

    // The new fluid density and velocity are obtained from the populations.
    momenta();

    return;
}

//****************************************Drag Force********************************************************
double drag(double den_c,double v_x){
  double C_D=0.45;
  return 0.5 * C_D * den_c * SQ(v_x) * 50.;
}

//*****************************FFT*************************************

void _fft(cplx buf[], cplx out[], int n, int step)
{
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);

		for (int i = 0; i < n; i += 2 * step) {
			cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}

void fft(cplx buf[], int n)
{
	cplx out[n];
	for (int i = 0; i < n; i++) out[i] = buf[i];

	_fft(buf, out, n, 1);
}

void show(const char *s, cplx buf[]) {
	printf("%s", s);
	for (int i = 0; i < Nx; i++)
		if (!cimag(buf[i]))
			printf("%g ", creal(buf[i]));
		else
			printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
}

//++++++++++++++++++++++++++++++init() & mask()++++++++++++++++++++++++++++++++++++++
//++++++init() inicializa el tipo de nodo de toda la red a 0+++++++++++++++++++++++++
//++++++mask() asigna el tipo de nodo (fluido=2, frontera=1, solido=0)+++++++++++++++

void init_type(){
  //alojar memoria para la densidad y las velocidades (vel_x y vel_y)
  type = new int*[Nx];

  for(int X = 0; X < Nx; X++) {
    type[X] = new int[Ny];
  }

  // Initialize the fluid density and velocity.
  for(int X = 0; X < Nx; X++) {
    for(int Y = 0; Y < Ny; Y++) {
      type[X][Y] = 0;
    }
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void mask(int x_0, int y_0, int r){
  for(int X=0; X<Nx; X++){
    for(int Y=0; Y<Ny; Y++){
      //printf("SQ(X-x_0)+SQ(Y-y_0)= %d\n", SQ((X-x_0))+SQ((Y-y_0)));
      if(SQ((X-x_0))+SQ((Y-y_0))<SQ(((double)r-epsilon))){
        type[X][Y]=0;
      }
      else if(SQ((X-x_0))+SQ((Y-y_0))>SQ(((double)r-epsilon))&& SQ((X-x_0))+SQ((Y-y_0))<SQ(r)) {
        type[X][Y]=1;
        fronteras++;
      }
      else type[X][Y]=2;
    }
  }
  return;
}

//++++++++++++++++++++++++sqared mask ++++++++++++++++++++++++++++++++++
  void mask_sq(int x_0, int y_0, int r){
    for(int X=0; X<Nx; X++){
      for(int Y=0; Y<Ny; Y++){
        type[X][Y]=2;
        if(X<(x_0+r+epsilon_sq) && X>(x_0-r-epsilon_sq)){
          if(Y<(y_0+r+epsilon_sq) && Y>(y_0-r-epsilon_sq)){
            type[X][Y]=1;
            fronteras++;
          }
        }
        if(X<(x_0+r-epsilon_sq) && X>(x_0-r+epsilon_sq)){
          if(Y<(y_0+r-epsilon_sq) && Y>(y_0-r+epsilon_sq)){
            type[X][Y]=0;
            fronteras--;
          }
        }
      }
    }
    return;
  }

//++++++++++++++++++++++++++++++++++++++++++++++++
void no_mask(){
  for(int X=0; X<Nx; X++){
    for(int Y=0; Y<Ny; Y++){
      type[X][Y]=2;
    }
  }
  return;
}

//+++++++++++++++++++++inv()+++++++++++++++++++++++++
//+++++Esta funciòn regresa la direcciòn inversa+++++

int inv(int j){
  int s=0;
  if(1<=j && j<=2){
    s=j+2;
  }
  else if(3<=j && j<=4){
    s=(((j+2))%4);
  }
  else if(4<j && j<7){
    s=j+2;
  }
  else{ s=j-2;}
  return s;
}

//++++++++++++++++++++++++++++++++++MEA()++++++++++++++++++++++++++++++++
//+++Calcula las fuerzas sobre fronteras sòlidas mediante el intercambio de momentum+++++
//++++dP=dx*dx Sum_(x_i^w) (f_i^in + f_(-i)^out) c_i+++++++++++++++++++++++++
//++++Aqui se utiliza para calcular las fuerzas de drag y lift sobre la barrera++++++++++
/*void MEA(){
  if(type[(X+cx[k])][(Y+cy[k])]==1){

  }
}*/
