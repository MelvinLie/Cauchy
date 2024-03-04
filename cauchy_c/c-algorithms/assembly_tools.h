#ifndef ASSEMBLYTOOLS_H
#define ASSEMBLYTOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss_legendre.h"

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif
#ifndef FABS
	#define FABS(a) ((a)>=0?(a):-(a))
#endif

/*
TO DO:
	- include segment length to discretization \Gamma, saves computation during the assembly.
	- adapt the scaling parameter R for unique solvability
	- include analytic integral evaluations, saves computation time and boosts accuracy for close elements.
	- compute multipoles directly by gaussian integration.

*/



/*--------------------------------------------------------------------
\Gamma is a discritization of the boundary into segments.
	o double *r:		Container for the N node coordinates of \Gamma.
						r = [x_0,y_0,x_1,y_1,...x_N,y_N,x_0,y_0]
						!!!It is important that the array closes the boundary
						That means that r ends with (x_0,y_0)!!!
	o int N:			Number of nodes. (without counting (x_0,y_0) twice)
	o int i:			Counting variable.
	o int j:			Counting variable.
	o double R:			Used to scale the argument in Green's function! 
						Important for unique solvability of the problem.
	o double *n:		Normal vectors on each segment.
						n = [n_x^0 , n_y^0 , ..., n_x^S , n_y^S]
	o double *r_eval:	Evaluation positions. Used for forward problems.
--------------------------------------------------------------------*/
typedef struct disc {

  double *r;
  int N;
  int i;
  int j;
  double R;
  double *n;
  double r_eval[2];

} disc;

/* Make a mesh parameterization*/
struct disc* make_parameterization(double* mesh, int N, double R);

/* Call gauss_legendre, increasing n until atol is reached*/
double integrate_1D(double (*f)(double,void*), void* data, double a, double b, double atol);

/* Call gauss_legendre, increasing n until atol is reached*/
double integrate_2D(double (*f)(double,double,void*), void* data, double a, double b, double c, double d, double atol);

/* Compute the lenght of a segment*/
double get_segment_lenght(struct disc *gamma, int seg);

/* Assemble the matrix corresponding to the single layer potential operator in a Galerkin scheme.*/
double*  assemble_V_c( disc *gamma);

/* Evaluate the kernel of the single layer potential operator.*/
double sl_kernel(double s, double t, void* data);

/* Evaluate the kernel of the double layer potential operator, multiplied with a rising linear basis function.*/
double dl_kernel_0(double s, double t, void* data);

/* Evaluate the kernel of the double layer potential operator, multiplied with a decreasing linear basis function.*/
double dl_kernel_1(double s, double t, void* data);

/* Compute the normal vectors on Gamma.*/
double* compute_normal_vectors(struct disc *gamma);

/* Assemble the matrix corresponding to the double layer potential operator in a Galerkin scheme.*/
double* assemble_K_c(struct disc *gamma);

/* Evaluate the gradient of the representation formula at positions given in r.*/
double* evaluate_grad_u(const double* r, const int M, struct disc *gamma, const double *v,const double *w);

/* Evaluate the x component of the kernen of the gradient of the single layer potential.*/
double eval_grad_sl_kernel_x(double s, void* data);

/* Evaluate the y component of the kernen of the gradient of the single layer potential.*/
double eval_grad_sl_kernel_y(double s, void* data);

/* Assemble the matrix M with [M]_{k,i} = <\varphi_k^0, \varphi_i^1 >_\Gamma, for constant \varphi_k^0 and linear \varphi_i^1*/
void assemble_M_c(struct disc *gamma, unsigned long long *rc, double *M);

double* assemble_D_c(struct disc *gamma,  double *V);

#endif
