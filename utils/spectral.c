#include "spectral.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

double
compute_two_norm(double **u_mat, int c_id, int size);
void
compute_project(double **a_mat, double **e_mat, 
		double *project, int c_id, 
		int diff, int size);
void 
minus(double **a_mat, double ** u_mat,
		double *project, int c_id,
		int size);
double
compute_inner(double **a_mat, double **e_mat,
		int size, int c_a, int c_e);
void
mat_mul(double **a_mat, double **q_mat,
		double **r_mat, int size);
double
compute_terminate(double **a_mat, int size);
void 
householder(double **a_mat, int size);
double 
compute_alpha(double **a_mat, int iter, int size);
void
compute_v(double **a_mat, double *v, 
		double alpha, double r, 
		int iter, int size);
void
multiply_and_minus(double *v, double **p, 
		double **I, int size);
void
print_mat(double **a_mat, int size);
int 
check_is_triangularized(double** a_mat, int size);
void 
init_mat(double **mat, int size);
void 
mat_copy(double **from, double **to, int size);
double
compute_gap(double *vector, int size);
void 
q_sort(double *numbers, int left, int right);
void 
compute_partition(double *vector, double alpha, 
			double beta, int *num_pos, 
			int *num_neg, int size, 
			double gap);

/*
 *return the gap
 **/
double
compute_second_eigen_vector(pdism d_mat, double *vector){	
	int i, j, k;
	int size = d_mat->num_spc;
	double *data = (double*)malloc(sizeof(double)*size*size);
	for(i=0; i<size; i++)
		for(j=0; j<size; j++)
			data[i*size+j] = d_mat->q_mat[i][j];
	gsl_matrix_view m 
		= gsl_matrix_view_array (data, size, size);

	gsl_vector *eval = gsl_vector_alloc (size);
	gsl_matrix *evec = gsl_matrix_alloc (size, size);

	gsl_eigen_symmv_workspace * w = 
		gsl_eigen_symmv_alloc (size);

	gsl_eigen_symmv (&m.matrix, eval, evec, w);

	gsl_eigen_symmv_free (w);

	gsl_eigen_symmv_sort (eval, evec, 
			GSL_EIGEN_SORT_ABS_ASC);

	//for (i = 0; i < size; i++)
	//{
	//	double eval_i 
	//		= gsl_vector_get (eval, i);
	//	gsl_vector_view evec_i 
	//		= gsl_matrix_column (evec, i);

	//	printf ("eigenvalue = %g\n", eval_i);
	//	printf ("eigenvector = \n");
	//	gsl_vector_fprintf (stdout, 
	//			&evec_i.vector, "%g");
	//}

	double max=-1000000, second_max=-1000000;
	int second =0;
	for (i = 0; i < size; i++){
		double eval_i = gsl_vector_get (eval, i);
		if(eval_i>max)
			max=eval_i;
	}
	for (i = 0; i < size; i++){
		double eval_i = gsl_vector_get (eval, i);
		if(eval_i>second_max && eval_i < max){
			second_max = eval_i;
			second = i;	
		}
	}
	gsl_vector_view evec_i = gsl_matrix_column (evec, second);
	//printf("second largest is %d %f\n", second, second_max);
	for(i=0; i<size; i++){
		vector[i] = gsl_matrix_get(evec, i, second);
	}

	gsl_vector_free (eval);
	gsl_matrix_free (evec);
	double gap = compute_gap(vector, size);
	q_sort(vector, 0, size-1);
	return gap;
	//double **e_mat = (double**)malloc(sizeof(double*)*d_mat->num_spc);
	//double **u_mat = (double**)malloc(sizeof(double*)*d_mat->num_spc);
	//double **q_mat = (double**)malloc(sizeof(double*)*d_mat->num_spc);
	//double **r_mat = (double**)malloc(sizeof(double*)*d_mat->num_spc);
	//double **p_q = (double**)malloc(sizeof(double*)*d_mat->num_spc);
	//double **tmp_pq = (double**)malloc(sizeof(double*)*d_mat->num_spc);
	//for(i=0; i<d_mat->num_spc; i++){
	//	e_mat[i] = (double*)malloc(sizeof(double)*d_mat->num_spc);
	//	u_mat[i] = (double*)malloc(sizeof(double)*d_mat->num_spc);
	//	q_mat[i] = (double*)malloc(sizeof(double)*d_mat->num_spc);
	//	r_mat[i] = (double*)malloc(sizeof(double)*d_mat->num_spc);
	//	p_q[i] = (double*)malloc(sizeof(double)*d_mat->num_spc);
	//	tmp_pq[i] = (double*)malloc(sizeof(double)*d_mat->num_spc);
	//	p_q[i][i] = 1.0;
	//}
	//double *project = (double*)malloc(sizeof(double)*d_mat->num_spc);
	//int terminate = FALSE;
	////perform triagulalization
	//print_mat(d_mat->q_mat, d_mat->num_spc);
	////exit(1);
	//householder(d_mat->q_mat, d_mat->num_spc);
	////perform QR algorithm by using Gram–Schmidt process
	//int count=0;
	//while(terminate == FALSE && count<1000){
	//	count++;
	//	init_mat(e_mat, d_mat->num_spc);
	//	init_mat(u_mat, d_mat->num_spc);
	//	init_mat(q_mat, d_mat->num_spc);
	//	init_mat(r_mat, d_mat->num_spc);
	//	//compute e_0 and u_0
	//	for(i=0; i<d_mat->num_spc; i++)
	//		u_mat[i][0] = d_mat->q_mat[i][0];
	//	double two_norm = compute_two_norm(u_mat, 0, d_mat->num_spc);
	//	for(i=0; i<d_mat->num_spc; i++)
	//		e_mat[i][0] = u_mat[i][0]/two_norm;	
	//	//for(i=0; i<d_mat->num_spc; i++){
	//	//	printf("%3.1f\t", e_mat[i][0]);
	//	//}
	//	//printf("\n");

	//	for(i=1; i<d_mat->num_spc; i++){			
	//		for(j=0; j<d_mat->num_spc; j++)
	//			project[j] = 0;
	//		for(j=0;j<i; j++)
	//			compute_project(d_mat->q_mat, e_mat, project, i, j+1, d_mat->num_spc);
	//		minus(d_mat->q_mat, u_mat, project, i, d_mat->num_spc);
	//		two_norm = compute_two_norm(u_mat, i, d_mat->num_spc);
	//		for(j=0; j<d_mat->num_spc; j++)
	//			e_mat[j][i] = u_mat[j][i]/two_norm;
	//		//printf("%3.1f ", two_norm);
	//	}
	//	//
	//	mat_copy(e_mat, q_mat, d_mat->num_spc);
	//	for(i=0; i<d_mat->num_spc; i++)
	//		for(j=i; j<d_mat->num_spc; j++)
	//			r_mat[i][j] = compute_inner(d_mat->q_mat, e_mat, d_mat->num_spc, j, i);
	//	mat_mul(d_mat->q_mat, q_mat, r_mat, d_mat->num_spc);
	//	mat_copy(p_q, tmp_pq, d_mat->num_spc);
	//	mat_mul(p_q, q_mat, tmp_pq, d_mat->num_spc);
	//	//print_mat(d_mat->q_mat, d_mat->num_spc);
	//	//print_mat(q_mat, d_mat->num_spc);
	//	//print_mat(r_mat, d_mat->num_spc);
	//	terminate = compute_terminate(d_mat->q_mat, d_mat->num_spc);
	//}
	////find the second largest eigen value
	//int second_idx = second_largest(d_mat->q_mat, d_mat->num_spc);
	//for(i=0; i<d_mat->num_spc; i++)
	//	vector[i] = r_mat[i][second_idx];
	//print_mat(d_mat->q_mat, d_mat->num_spc);
	//print_mat(p_q, d_mat->num_spc);
	//exit(1);
	////frees
	//for(i=0; i<d_mat->num_spc; i++){
	//	free(e_mat[i]);
	//	free(u_mat[i]);
	//	free(q_mat[i]);
	//	free(r_mat[i]);
	//	free(p_q[i]);
	//	free(tmp_pq[i]);
	//}
	//free(e_mat);
	//free(u_mat);
	//free(q_mat);
	//free(r_mat);
	//free(project);
	//free(p_q);
	//free(tmp_pq);
}

double
compute_gap(double *vector, int size){
	int i, j, k;
	double min_pos = 10000000;
	double max_neg = -10000000;
	for(i=0; i<size; i++)
		if(vector[i]<min_pos && vector[i]>0)
			min_pos = vector[i];
	for(i=0; i<size; i++)
		if(vector[i]>max_neg && vector[i]<=0)
			max_neg = vector[i];
	return (min_pos - max_neg);
}

void 
mat_copy(double **from, double **to, int size){
	int i, j;
	for(i=0; i<size; i++)
		for(j=0; j<size; j++)
			to[i][j] = from[i][j];
}

void 
init_mat(double **mat, int size){
	int i, j;
	for(i=0; i<size; i++)
		for(j=0; j<size; j++)
			mat[i][j] = 0;
}

int 
second_largest(double **a_mat, int size){
	int i, j, k;
	double max=0;
	for(i=0; i<size; i++)
		if(a_mat[i][i]>max)
			max = a_mat[i][i];
	double second_max =0;
	int second_idx =1;
	for(i=0; i<size; i++)
		if(a_mat[i][i]>second_max && a_mat[i][i]<max){
			second_max = a_mat[i][i];
			second_idx = i;
		}
	return second_idx;	
}

double
compute_two_norm(double **u_mat, int c_id, int size){
	int i, j, k;
	double sum =0;
	for(i=0; i<size; i++){
		sum += u_mat[i][c_id]*u_mat[i][c_id];
	}
	return sqrt(sum);
}

void
compute_project(double **a_mat, double **e_mat, 
		double *project, int c_id, 
		int diff, int size){
	int i, j, k;
	//compute the inner product
	double inner_ea = 0;
	double inner_ee = 0;
	for(i=0; i<size; i++)
		inner_ea += a_mat[i][c_id] * e_mat[i][c_id-diff];
	//for(i=0; i<size; i++)
	//	inner_ee += e_mat[i][c_id-diff]*e_mat[i][c_id-diff];
	//for(i=0; i<size; i++)
	//	project[i] += (e_mat[i][c_id-diff]*inner_ea)/inner_ee;
	for(i=0; i<size; i++)
		project[i] += (e_mat[i][c_id-diff]*inner_ea);
}

void 
minus(double **a_mat, double ** u_mat,
		double *project, int c_id,
		int size){
	int i;
	for(i=0; i<size; i++)
		u_mat[i][c_id] = a_mat[i][c_id]-project[i];
}

double
compute_inner(double **a_mat, double **e_mat,
		int size, int c_a, int c_e){
	int i;
	double sum = 0;
	for(i=0; i<size; i++){
		sum += a_mat[i][c_a]*e_mat[i][c_e];
		//printf("a:%f e:%f sum:%f \n", a_mat[i][c_a], e_mat[i][c_e], sum);
	}
	return sum;
}

void
mat_mul(double **a_mat, double **q_mat,
		double **r_mat, int size){
	int i, j, k;
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			a_mat[i][j]=0;
			for(k=0;k<size;k++){
				a_mat[i][j]=a_mat[i][j]+(r_mat[i][k] * q_mat[k][j]);
			}
		}
	}
}

double
compute_terminate(double **a_mat, int size){
	int i, j;
	int result = TRUE;
	double thresh = 0.00001;
	for(i=0; i<size; i++)
		for(j=i+1; j<size; j++){
			//printf("i %d j %d val %f\n", i, j, fabs(a_mat[j][i]));
			if(fabs(a_mat[j][i])>thresh){
				result = FALSE;
				//printf("i %d j %d val %f\n", i, j, a_mat[j][i]);
			}
		}
	//printf("\n");
	return result;
}

void householder(double **a_mat, int size){
	int i, j, k;
	//p is the household matrix
	double **p = (double**)malloc(sizeof(double*)*size);
	double **I = (double**)malloc(sizeof(double*)*size);
	double **tmp = (double**)malloc(sizeof(double*)*size);
	double *v = (double*)malloc(sizeof(double)*size);
	for(i=0; i<size; i++){
		p[i] = (double*)malloc(sizeof(double)*size);
		I[i] = (double*)malloc(sizeof(double)*size);
		tmp[i] = (double*)malloc(sizeof(double)*size);
		I[i][i] = 1;
	}
	for(i=0; i<size-1; i++){
		int sign = a_mat[i+1][i]>0?-1:1;
		double alpha = sign*compute_alpha(a_mat, i, size);
		double r = sqrt((alpha*alpha - a_mat[i+1][i]*alpha)/2);
		//printf("sign: %d alpha %3.2f r %3.2f\n", sign, alpha, r);
		compute_v(a_mat, v, alpha, r, i, size);
		multiply_and_minus(v, p, I, size);
		mat_mul(tmp, p, a_mat, size);
		mat_mul(a_mat, tmp, p, size);
		if(check_is_triangularized(a_mat, size)==TRUE)
			break;
	}
	//free
	for(i=0; i<size; i++){
		free(p[i]);
		free(I[i]);
		free(tmp[i]);
	}
	free(p);
	free(I);
	free(tmp);
	free(v);
}

int 
check_is_triangularized(double** a_mat, int size){
	int i, j, k;
	for(i=0; i<size; i++){
		int left_pos = i>1?i-2:-1;
		int right_pos = i<size-2?i+2:size;
		for(j=0; j<left_pos; j++)
			if(fabs(a_mat[i][j]) > 0.000001){
				return FALSE;
			}
		for(j=right_pos; j<size; j++)
			if(fabs(a_mat[i][j]) > 0.000001){
				return FALSE;
			}
	}
	return TRUE;
}

void
print_mat(double **a_mat, int size){
	int i, j, k;
	for(i=0; i<size; i++){
		for(j=0; j<size; j++){
			printf("%4.3f\t", a_mat[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
}

double 
compute_alpha(double **a_mat, int iter, int size){
	int i;
	double sum =0;
	for(i=iter+1; i<size; i++){
		sum += a_mat[i][iter]*a_mat[i][iter];
	}
	return sqrt(sum);
}

void
compute_v(double **a_mat, double *v, 
		double alpha, double r, 
		int iter, int size){
	int i,j,k;
	for(i=0; i<=iter; i++)
		v[i] = 0;
	v[iter+1] = (a_mat[iter+1][iter]-alpha)/(2*r);
	for(i=iter+2; i<size; i++)
		v[i] = a_mat[i][iter]/(2*r);
}

void
multiply_and_minus(double *v, double **p, double **I, int size){
	int i, j;
	for(i=0; i<size; i++)
		for(j=0; j<size; j++)
			p[i][j] = 2*v[i]*v[j];
	for(i=0; i<size; i++)
		for(j=0; j<size; j++)
			p[i][j] = I[i][j]-p[i][j];

}

void 
q_sort(double *numbers, int left, int right)
{
        double pivot;
	int l_hold, r_hold;

        l_hold = left;
        r_hold = right;
        pivot = numbers[left];
        while (left < right)
        {
                while ((numbers[right] >= pivot) && (left < right))
                        right--;
                if (left != right)
                {
                        numbers[left] = numbers[right];
                        left++;
                }
                while ((numbers[left] <= pivot) && (left < right))
                        left++;
                if (left != right)
                {
                        numbers[right] = numbers[left];
                        right--;
                }
        }
        numbers[left] = pivot;
        pivot = left;
        left = l_hold;
        right = r_hold;
        if (left < pivot)
                q_sort(numbers, left, pivot-1);
        if (right > pivot)
                q_sort(numbers, pivot+1, right);
}

void 
compute_partition(double *vector, double alpha, 
			double beta, int *num_pos, 
			int *num_neg, int size, 
			double gap){
	int i, j, k;
	if(gap > alpha){
		for(i=0; i<size; i++)
			if(vector[i]< 0)
				(*num_neg) ++;
		*num_pos = size - *num_neg;
		if((*num_pos) <= 1)
			*num_pos = 2;
		if((*num_neg) <= 1)
			*num_neg = 2;
	}
	else{
		int *is_overlap = (int*)malloc(sizeof(int)*size);
		for(i=0; i<size; i++)
			is_overlap[i] = FALSE;
		int bound_pos = -1, bound_neg = -1;
		for(i=0; i<size-1; i++)
			if(vector[i]<0 && vector[i+1]>=0)
				bound_neg = i;
		for(i=size-1; i>0; i--)
			if(vector[i]>=0 && vector[i-1]<0)
				bound_pos = i;
		for(i=bound_neg; i>0; i--){
			if((vector[i]-vector[i-1]) < beta)
				is_overlap = TRUE;
			else
				break;
		}
		for(i=bound_pos; i<size-1; i++){
			if((vector[i+1]-vector[i]) < beta)
				is_overlap = TRUE;
			else
				break;
		}
		for(i=0; i<size; i++)
			if(vector[i]<0 || is_overlap[i]==TRUE)
				*num_neg++;
		for(i=0; i<size; i++)
			if(vector[i]>=0 || is_overlap[i]==TRUE)
				*num_pos++;
		free(is_overlap);
	}
}
