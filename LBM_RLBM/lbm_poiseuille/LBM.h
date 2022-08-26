/* This code accompanies
 *   The Lattice Boltzmann Method: Principles and Practice
 *   T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
 *   ISBN 978-3-319-44649-3 (Electronic)
 *        978-3-319-44647-9 (Print)
 *   http://www.springer.com/978-3-319-44647-9
 *
 * This code is provided under the MIT license. See LICENSE.txt.
 *
 * Author: Orest Shardt
 *
 */
#ifndef __LBM_H
#define __LBM_H

const unsigned int scale = 1;
const unsigned int NX = 5*scale;     // Dimensiones de la red
const unsigned int NY = NX;

const unsigned int ndir = 9;
const size_t mem_size_0dir   = sizeof(double)*NX*NY;           // calcula el tamaño
const size_t mem_size_n0dir  = sizeof(double)*NX*NY*(ndir-1);  // en bytes de las n direcciones
const size_t mem_size_scalar = sizeof(double)*NX*NY;           // y de un escalar

const double w0 = 4.0/9.0;  // zero weight
const double ws = 1.0/9.0;  // adjacent weight
const double wd = 1.0/36.0; // diagonal weight

const double tau=sqrt(3/16)+0.5; // en matlab (Chapter 5, tau=sqrt(3/16)+0.5 gives exact solution)  relaxation time (BGK model)
const double nu=(2*tau-1)/6;         // kinematic shear viscosity
const double Re=NY*u_max/nu;
const double omega = 1./nu;    // omega=1/tnu en Chapter 5

// Taylor-Green parameters
const double u_max = 0.1/scale;   // velocidad maxima en la red
const double rho0 = 1.0;           // densidad inicial

const unsigned int NSTEPS = 10000*scale*scale;    // Numero de pasos que se realizarà el metodo
const unsigned int NSAVE  =  50*scale*scale;    // Numeros de referencia para el guardado de datos de salida
const unsigned int NMSG   =  50*scale*scale;

// compute L2 error and energy?
// disable for speed testing
const bool computeFlowProperties = true;

// suppress verbose output
const bool quiet = true;

//Aqui se declaran las funciones que mas adelante se definen (o en otros archivos)
void taylor_green(unsigned int,unsigned int,unsigned int,double*,double*,double*);
void taylor_green(unsigned int,double*,double*,double*);
void stream_collide_save(double*,double*,double*,double*,double*,double*,bool);
void init_equilibrium(double*,double*,double*,double*,double*);
void compute_flow_properties(unsigned int,double*,double*,double*,double*);
void report_flow_properties(unsigned int,double*,double*,double*);
void save_scalar(const char*,double*,unsigned int);

// Funciòn que calcula la posiciòn lineal de las poblaciones con velocidad c_0
inline size_t field0_index(unsigned int x, unsigned int y)
{
    return NX*y+x;
}

// Funciòn que calcula la posiciòn lineal de los escalares en las poblaciones
inline size_t scalar_index(unsigned int x, unsigned int y)
{
    return NX*y+x;
}

// Funciòn que calcula la posiciòn lineal de las poblaciones con velocidad c_{1-8}
inline size_t fieldn_index(unsigned int x, unsigned int y, unsigned int d)
{
    return (ndir-1)*(NX*y+x)+(d-1);
}

#endif /* __LBM_H */
