cdef extern from "assembly_tools.h":

	ctypedef struct disc:
		pass

	disc* make_parameterization(double* mesh, int N, double R)

	double integrate_1D(double (*f)(double, void*), void* data, double a, double b, double atol)

	double integrate_2D(double (*f)(double, double, void*), void* data, double a, double b, double c, double d, double atol)

	double get_segment_lenght(disc *gamma, int seg)

	double*  assemble_V_c( disc *gamma)

	double sl_kernel(double s, double t, void* data)

	double dl_kernel_0(double s, double t, void* data)

	double dl_kernel_1(double s, double t, void* data)

	double* compute_normal_vectors(disc *gamma)

	double* assemble_K_c(disc *gamma)

	double* evaluate_grad_u(const double* r, const int M, disc *gamma, const double *v,const double *w)

	double eval_grad_sl_kernel_x(double s, void* data)

	double eval_grad_sl_kernel_y(double s, void* data)

	void assemble_M_c(disc *gamma, unsigned long long *rc, double *M)

	double* assemble_D_c(disc *gamma,  double *V)