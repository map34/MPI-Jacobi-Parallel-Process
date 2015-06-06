/*
 * Mochamad Prananda
 * HW7
 */



#include "MyHeader.h"



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CreateMatrix
 *  Description:  
 * =====================================================================================
 */
SparseMatrix* CreateMatrix (int my_rank, int n, int comm_sz)
{
	return NULL;
}		/* -----  end of function CreateMatrix  ----- */


/* ---- Sparse Matrix Multiplication --- */
void MatVec(const struct SparseMatrix *local_A, double *local_x, double *local_y, int n, int comm_sz, MPI_Comm comm)
{

	int local_i, k,i;
	int local_n = n/comm_sz;

	double *globalX = (double*)malloc(sizeof(double)*n);

	

	MPI_Allgather(local_x, local_n, MPI_DOUBLE, globalX, local_n, MPI_DOUBLE, comm);

/*for (i = 0; i < n; i++)
	{
		printf("x=%e\n", globalX[i]); 
	}*/

	for (local_i = 0; local_i < local_n; ++local_i)
	{
		double sum = 0.0;
		for (k = local_A->iRow[local_i]; k < local_A->iRow[local_i+1]; ++k)
		{
			sum += local_A->values[k] * globalX[local_A->jCol[k]];
		}
		local_y[local_i] = sum;
	}
	free(globalX);

}


/* ---- Display Sparse Matrix --- */
void Printf(FILE *stream, const struct SparseMatrix *A, int comm_sz)
{

	int i, k;
	for (i = 0; i <A->m/comm_sz ; ++i)
	{
		for (k = A->iRow[i]; k < A->iRow[i+1]; ++k)
		{
			fprintf(stream, "%d %d %e\n", i, A->jCol[k], A->values[k]);
			
		}
	}

}


/* --- Output a vector --- */
void PrintVector(FILE *stream, const struct Vector *b)
{

	int i;
	for (i = 0; i < b->n; ++i)
	{
		fprintf(stream, "%e\n", b->values[i]);
	}

}


/* --- Function to get current time --- */
double getTime() 
{

	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec/1000000.0;

}


