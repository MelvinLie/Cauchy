#include "assembly_tools.h"

struct disc* make_parameterization(double* mesh, int N, double R){

	// allocate the memory for a discretization
	disc *gamma	= (disc*) calloc(1, sizeof(disc));

	// fill it
	gamma->r = mesh;
	gamma->N = N;
	gamma->i = 0;
	gamma->j = 0;
	gamma->R = R;
	gamma->n = compute_normal_vectors(gamma);

	return gamma;

}

double eval_grad_sl_kernel_x(double s, void* data){

	disc *gamma = (disc*) data;

	double d_x = gamma->r_eval[0] - s*(gamma->r[2*gamma->i+2] - gamma->r[2*gamma->i])   -  gamma->r[2*gamma->i];
	double d_y = gamma->r_eval[1] - s*(gamma->r[2*gamma->i+3] - gamma->r[2*gamma->i+1]) -  gamma->r[2*gamma->i+1];

	return d_x/(d_x*d_x + d_y*d_y);

}

double eval_grad_sl_kernel_y(double s, void* data){

	disc *gamma = (disc*) data;

	double d_x = gamma->r_eval[0] - s*(gamma->r[2*gamma->i+2] - gamma->r[2*gamma->i])   -  gamma->r[2*gamma->i];
	double d_y = gamma->r_eval[1] - s*(gamma->r[2*gamma->i+3] - gamma->r[2*gamma->i+1]) -  gamma->r[2*gamma->i+1];

	return d_y/(d_x*d_x + d_y*d_y);

}

double sl_kernel(double s, double t, void* data)
{

	disc *gamma = (disc*) data;


	double num_x = t*(gamma->r[2*gamma->j+2] - gamma->r[2*gamma->j]) + gamma->r[2*gamma->j]
								 - s*(gamma->r[2*gamma->i+2] - gamma->r[2*gamma->i]) -  gamma->r[2*gamma->i];

	double num_y = t*(gamma->r[2*gamma->j+3] - gamma->r[2*gamma->j+1]) + gamma->r[2*gamma->j+1]
								 - s*(gamma->r[2*gamma->i+3] - gamma->r[2*gamma->i+1]) -  gamma->r[2*gamma->i+1];

	//printf("%f\n",log(sqrt(num_x*num_x + num_y*num_y)/gamma->R));
	return log(sqrt(num_x*num_x + num_y*num_y)/gamma->R);
	//return 0.;
}

double dl_kernel_0(double s, double t, void* data)
{

	struct disc *gamma = (struct disc*) data;
	double dx, dy;

	dx = t*(gamma->r[2*gamma->i+2] - gamma->r[2*gamma->i]) + gamma->r[2*gamma->i]
     - s*(gamma->r[2*gamma->j+2] - gamma->r[2*gamma->j]) -  gamma->r[2*gamma->j];

	dy = t*(gamma->r[2*gamma->i+3] - gamma->r[2*gamma->i+1]) +  gamma->r[2*gamma->i+1]
	   - s*(gamma->r[2*gamma->j+3] - gamma->r[2*gamma->j+1]) -  gamma->r[2*gamma->j+1];

	//printf("%f\n",1/(dx*dx + dy*dy));
	return s*(gamma->n[2*gamma->j]*dx+ gamma->n[2*gamma->j+1]*dy)/(dx*dx + dy*dy);
	//return 0.;
}

double dl_kernel_1(double s, double t, void* data)
{

	struct disc *gamma = (struct disc*) data;
	double dx, dy;

	dx = t*(gamma->r[2*gamma->i+2] - gamma->r[2*gamma->i]) + gamma->r[2*gamma->i]
     - s*(gamma->r[2*gamma->j+2] - gamma->r[2*gamma->j]) -  gamma->r[2*gamma->j];

	dy = t*(gamma->r[2*gamma->i+3] - gamma->r[2*gamma->i+1]) +  gamma->r[2*gamma->i+1]
	   - s*(gamma->r[2*gamma->j+3] - gamma->r[2*gamma->j+1]) -  gamma->r[2*gamma->j+1];

	//printf("%f\n",1/(dx*dx + dy*dy));
	//printf("%f\n",log(sqrt(num_x*num_x + num_y*num_y)/gamma->R));
	return (1-s)*(gamma->n[2*gamma->j]*dx+ gamma->n[2*gamma->j+1]*dy)/(dx*dx + dy*dy);
	//return 0.;
}

double get_segment_lenght(struct disc *gamma, int seg){

	int seg_indx = seg*2;
	double len;

	if (seg < gamma->N-1){
		len = sqrt((gamma->r[seg_indx+2] - gamma->r[seg_indx])*(gamma->r[seg_indx+2] - gamma->r[seg_indx])	+ (gamma->r[seg_indx+3] - gamma->r[seg_indx+1])*(gamma->r[seg_indx+3] - gamma->r[seg_indx+1]));
	}
	else{
		len = sqrt((gamma->r[0] - gamma->r[seg_indx])*(gamma->r[0] - gamma->r[seg_indx])	+ (gamma->r[1] - gamma->r[seg_indx+1])*(gamma->r[1] - gamma->r[seg_indx+1]));
	}
	return len;

}

double integrate_1D(double (*f)(double,void*), void* data, double a, double b, double atol)
{

	/* error between two iterations*/
	double error;
	double approx_pre, approx_post;

	/*iteration variable*/
  int i = 3;

	/*convergence flag*/
	int convergence = 0;

	approx_pre = gauss_legendre(2,f,data,a,b);

	while (convergence == 0)
	{
		approx_post = gauss_legendre(i, f, data, a, b);
		//printf("i = %d, approx = %.15f\n",i,approx_post);
		error = approx_post-approx_pre;

		if(FABS(error) < atol)	
			convergence = 1;
		if(i > 1024)	{

			convergence = 1;
			// printf("Gaussian Integration: Maximum iterations reached!\n");

		}

		approx_pre = approx_post;
		i = i + 1;

	}

	return approx_post;

}

double integrate_2D(double (*f)(double,double,void*), void* data, double a, double b, double c, double d, double atol){

	/* error between two iterations*/
	double error;
	double approx_pre, approx_post;

	/*iteration variable*/
	int i = 3;

	/*convergence flag*/
	int convergence = 0;

	approx_pre = gauss_legendre_2D_cube(2,f,data,a,b,c,d);

	while (convergence == 0)
	{
		approx_post = gauss_legendre_2D_cube(i,f,data,a,b,c,d);
		//printf("iter = %d, approx = %.15f\n",i,approx_post);
		error = approx_post-approx_pre;

		if(FABS(error) < atol)		convergence = 1;
		if(i > 128)	{

			convergence = 1;
			// printf("Warning: Maximum iterations reached!\n");

		}

		approx_pre = approx_post;
		i = i + 1;

	}

	return approx_post;

}

double* test(){

	double *K;

	printf("Hello\n");
	K = (double*) malloc(sizeof(double));
	return K;

}


double* assemble_V_c(struct disc* gamma){

  int i, j;
  double *V,li,lj;

	V = (double*) calloc(gamma->N * gamma->N, sizeof(double)); // We could exploit symmetry already here. We could reduce the memory to 0.5 * N * ( N + 1 ) entries

	for(i = 0 ; i < gamma->N ; i++){
		for(j = i ; j < gamma->N ; j++){

			if (i == j){

				li = get_segment_lenght(gamma,i);

				V[i*gamma->N + j] = -li*li/2./PI*(log(li)-3./2.);

			}
			else{

				li = get_segment_lenght(gamma,i);
				lj = get_segment_lenght(gamma,j);	//CHECK!!!

				gamma->i = i;
				gamma->j = j;

				V[i*gamma->N + j] = -li*lj/2./PI*integrate_2D(sl_kernel, gamma, 0., 1., 0., 1.,1e-10);
				V[j*gamma->N + i] = V[i*gamma->N + j];

			}
		}
	}

	gamma->i = 0;
	gamma->j = 0;
  return V;

}


double* compute_normal_vectors(struct disc *gamma){

  int k;
	double *n;
	double norm;
	int K = gamma->N;

	n = (double*) calloc( K*2 ,sizeof(double));

  for(k = 0 ; k < K ; k++){

		if (k == K-1){

			n[2*K-2] = gamma->r[1] - gamma->r[2*K-1];		//nx = y_0 - y_{N-1}
			n[2*K-1] = gamma->r[2*K-2] - gamma->r[0];		//ny = x_{N-1} - x_0
			//norm
			norm = sqrt(n[2*K-2]*n[2*K-2] + n[2*K-1]*n[2*K-1]);
			n[2*K-2] *= 1/norm;
			n[2*K-1] *= 1/norm;

		}
		else{

			n[2*k] = gamma->r[2*(k+1)+1] - gamma->r[2*k+1];	//nx = y_{k+1} - y_k
			n[2*k+1] = gamma->r[2*k] - gamma->r[2*(k+1)];		//ny = x_k - x_{k+1}
			//norm
			norm = sqrt(n[2*k]*n[2*k] + n[2*k+1]*n[2*k+1]);
			n[2*k] *= 1/norm;
			n[2*k+1] *= 1/norm;

		}


	}
	return n;

}





double* assemble_K_c(struct disc *gamma){

  int i, j, j_pre, j_post;
  double *K,li,lj;




	K = (double*) calloc(gamma->N * gamma->N,sizeof(double));


	for(i = 0 ; i < gamma->N ; i++){
		for(j = 0 ; j < gamma->N ; j++){

			//printf("i = %d , j = %d\n",i,j);

			li = get_segment_lenght(gamma,i);

			if (j == 0){
				j_pre = gamma->N-1;
				j_post = 0;
			}
			else{
				j_pre = j-1;
				j_post = j;
			}

			//printf("j_pre = %d , j_post = %d\n",j_pre,j_post);

			if (i != j_pre){

				gamma->i = i;
				gamma->j = j_pre;

				lj = get_segment_lenght(gamma,j_pre);
				//integrate dl_kernel
				//printf("i = %d, j = %d, I = %f\n",i,j_pre,integrate_2D(dl_kernel_0, gamma, 0., 1., 0., 1.,1e-10));
				K[i*gamma->N + j] += li*lj/2./PI*integrate_2D(dl_kernel_0, gamma, 0., 1., 0., 1.,1e-10);
				//printf("I = %.12f\n",li*lj/2./PI*integrate_2D(dl_kernel_0, gamma, 0., 1., 0., 1.,1e-10));

			}
			if (i != j_post){

				gamma->i = i;
				gamma->j = j_post;

				lj = get_segment_lenght(gamma,j_post);
				//integrate dl dl_kernel
				//printf("i = %d, j = %d, I = %f\n",i,j_post,integrate_2D(dl_kernel_1, gamma, 0., 1., 0., 1.,1e-10));
				K[i*gamma->N + j] += li*lj/2./PI*integrate_2D(dl_kernel_1, gamma, 0., 1., 0., 1.,1e-10);
				//printf("I = %.12f\n",li*lj/2./PI*integrate_2D(dl_kernel_1, gamma, 0., 1., 0., 1.,1e-10));

			}
			//printf("V[%d,%d] = %.12f\n",i,j,V[i*gamma->N + j]);
		}
	}

  return K;

}

double* evaluate_grad_u(const double* r, const int M, struct disc *gamma, const double *v, const double *w){

	int k,m;
	double *grad;
	double d_sl_x, d_sl_y;
	double l_k, l_k_m1, curl_v;
	int k_m1, k_p1;

	grad = (double*) calloc(2*M,sizeof(double));



	for(m = 0; m < M ; m++){

		gamma->r_eval[0] = r[2*m];
		gamma->r_eval[1] = r[2*m+1];

		for(k = 0 ; k < gamma->N ; k++){

			gamma->i = k;
			l_k = get_segment_lenght(gamma,k);

			if (k == 0){
				l_k_m1 = get_segment_lenght(gamma,gamma->N-1);
				k_m1 = gamma->N-1;
			}
			else{
				l_k_m1 = get_segment_lenght(gamma,k-1);
				k_m1 = k-1;
			}

			if (k == gamma->N-1){
				k_p1 = 0;
			}
			else{
				k_p1 = k+1;
			}

			//TO DO: write Gaussian integration for double pointer to spare two integrations!
			//Integrate gradient of single layer potential
			d_sl_x = integrate_1D(eval_grad_sl_kernel_x, gamma, 0, 1, 1e-10);
			d_sl_y = integrate_1D(eval_grad_sl_kernel_y, gamma, 0, 1, 1e-10);


			grad[2*m] += l_k*w[k]*d_sl_x;
		 	grad[2*m+1] += l_k*w[k]*d_sl_y;


			//Integrate gradient of double layer potential we make use of integration by parts
			//curl_v = gamma->n[2*k_m1]*(gamma->r[2*k+1] - gamma->r[2*k_m1+1]) - gamma->n[2*k_m1+1]*(gamma->r[2*k] - gamma->r[2*k_m1]);
			curl_v = gamma->n[2*k]*(gamma->r[2*k+1] - gamma->r[2*k_p1+1]) - gamma->n[2*k+1]*(gamma->r[2*k] - gamma->r[2*k_p1]);


			grad[2*m] += (v[k_p1]-v[k]) * curl_v * d_sl_y / l_k;
			grad[2*m+1] -= (v[k_p1]-v[k]) * curl_v * d_sl_x / l_k;

			//curl_v = gamma->n[2*k]*(gamma->r[2*k+1] - gamma->r[2*k_p1+1]) - gamma->n[2*k+1]*(gamma->r[2*k] - gamma->r[2*k_p1]);

			//grad[0] += v[k] * curl_v * d_sl_y / l_k;
			//grad[1] -= v[k] * curl_v * d_sl_x / l_k;

		}

		grad[2*m] *= -1/2./PI;
		grad[2*m+1] *= -1/2./PI;

	}

	return grad;

}


void assemble_M_c(struct disc *gamma, unsigned long long *rc, double *M) {

	int k, i;



	for(k = 0 ; k < gamma->N; k++){

		M[2*k] = 0.5*get_segment_lenght(gamma,k);
        rc[4*k] = k;
		rc[4*k+1] = k;


		M[2*k+1] = M[2*k];
		rc[4*k+2] = k;

		if (k < gamma->N - 1) {
			rc[4*k+3] = k+1;
		} else {
			rc[4*k+3] = 0;
		}
    }

}




double* assemble_D_c(struct disc *gamma, double *V){

	double ai,bi,aj,bj,l,*D;
	int i,j,i_pre,i_post,j_pre,j_post;

	D = (double*) calloc(gamma->N * gamma->N,sizeof(double));

	if(V == NULL){
		//We need to assemble V first
		V = (double*) calloc(gamma->N * gamma->N, sizeof(double));
		V = assemble_V_c(gamma);

	}

	for(i = 0 ; i < gamma->N ; i++){

		if(i == 0) i_pre = gamma->N-1;
		else i_pre = i-1;

		if(i == gamma->N-1) i_post = 0;
		else i_post = i+1;

		l = get_segment_lenght(gamma,i_pre);
		//printf("l_i-1 = %f\n",l);
		ai = (gamma->n[2*i_pre]*(gamma->r[2*i+1]-gamma->r[2*i_pre+1]) - gamma->n[2*i_pre+1]*(gamma->r[2*i]-gamma->r[2*i_pre])) / l / l;

		l = get_segment_lenght(gamma,i);
		//printf("l_i = %f\n",l);
		bi = (gamma->n[2*i]*(gamma->r[2*i+1]-gamma->r[2*i_post+1]) - gamma->n[2*i+1]*(gamma->r[2*i]-gamma->r[2*i_post])) / l / l;

		for(j = i ; j < gamma->N; j++){

			if(j == 0) j_pre = gamma->N-1;
			else j_pre = j-1;

			if(j == gamma->N-1) j_post = 0;
			else j_post = j+1;

			//TO DO We can optimize run time, we will compute l twice here between two iterations which is not needed!
			l = get_segment_lenght(gamma,j_pre);
			//printf("l_j-1 = %f\n",l);
			aj = (gamma->n[2*j_pre]*(gamma->r[2*j+1]-gamma->r[2*j_pre+1]) - gamma->n[2*j_pre+1]*(gamma->r[2*j]-gamma->r[2*j_pre])) / l / l;

			l = get_segment_lenght(gamma,j);
			//printf("l_j = %f\n",l);
			bj = (gamma->n[2*j]*(gamma->r[2*j+1]-gamma->r[2*j_post+1]) - gamma->n[2*j+1]*(gamma->r[2*j]-gamma->r[2*j_post])) / l / l;

			//printf("i_pre = %d\n",i_pre);
			//printf("i_post = %d\n",i_post);
			//printf("j_pre = %d\n",j_pre);
			//printf("i_post = %d\n",j_post);
			//printf("ai = %f\n",ai);
			//printf("bi = %f\n",bi);
			//printf("aj = %f\n",aj);
			//printf("bj = %f\n",bj);
			D[i*gamma->N+j] = ai*aj*V[i_pre*gamma->N+j_pre] + ai*bj*V[i_pre*gamma->N+j] + bi*aj*V[i*gamma->N+j_pre] + bi*bj*V[i*gamma->N+j];
			D[j*gamma->N+i] = D[i*gamma->N+j];

		}
	}
	return D;

}
