/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  Jacobi Iteration using MPI
 *
 *        Version:  1.0
 *        Created:  05/31/2015 20:54:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mochamad Prananda (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include "MyHeader.h"
int main(int argc, char** argv)	
{
	int pNum;
	int myRank;
	int dim;
	int numRow;
	int rowOrder;
	FILE* out;

	SparseMatrix* A_local;
	double* diagInv_local;
	double* x_global;
	double* b_local; 
	double* r_local;
	double* y_local;

	double tmpNormB;
	double normB_global;
	double tmpNormR;
	double normR_global;

	/* *** START CODING *** */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	// Checks number of arguments
	if (argc <= 1)
	{
		if (myRank == 0)
		{
			fprintf(stderr,"Error: needs an argument\n");
		}
		MPI_Finalize();
		return -1;
	}

	// Checks the correct type of input
	char* ptr;
	dim = (int) strtol(argv[1],	&ptr, 10);
	if (*ptr != '\0')
	{
		if (myRank == 0)
		{
			fprintf(stderr,"Error: needs an integer\n");
		}
		MPI_Finalize();
		return -1;
	}

	
	numRow = dim / pNum;
	rowOrder = myRank * numRow;
	// get matrix
	int i,j;
	A_local = CreateMatrix(myRank, dim, pNum);

	// Get x values ** GLOBAL **
	x_global = (double*) malloc (sizeof(double)*dim);
	for (i = 0; i < dim; i++)
	{
		x_global[i] = 0.0;
	}

	// Get b values ** LOCAL **
	b_local = (double*) malloc (sizeof(double)*numRow);
	for (i = 0; i < numRow; i++)
	{
		b_local[i] = 1.0;
	}
	
	// JACOBI MPIi
	


	// Get diagonal matrix ** LOCAL **
	diagInv_local = (double*) malloc (sizeof(double)*numRow);
	
	for (i = 0; i < numRow; i++)
	{
		diagInv_local[i] = 0.5;
	}

	// Get norm of b
	for (i = 0; i < numRow; i++)
	{
		normB_global += b_local[i] *b_local[i];
	}
	tmpNormB = normB_global;

	MPI_Allreduce(&tmpNormB, &normB_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// Initialize y and r
	r_local = (double*) malloc (sizeof(double)*numRow);
	y_local = (double*) malloc (sizeof(double)*numRow);

	// JACOBI Process
	int l;
	double sum_local;
	for (i = 1; i < 2000000; i++)
	{
		// y = A*x
		for (j = 0; j < numRow; j++)
		{
			sum_local = 0.0;
			for (l = A_local->iRow[j]; l < A_local->iRow[j+1]; l++)
			{
				sum_local += A_local->values[l] * x_global[A_local->jCol[l]];
			}
			y_local[j] = sum_local;
		}

		// r = b - y, y = A*x
		for (j = 0; j < numRow; j++)
		{
			r_local[j] = b_local[j] - y_local[j];
		}
		
		// norm R calculations
		normR_global = 0.0;
		for (j = 0; j < numRow; j++)
		{
			normR_global += r_local[j]*r_local[j];
		}
		tmpNormR = normR_global;

		MPI_Allreduce(&tmpNormR, &normR_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// Break off point
		if (normR_global < 1.0e-8 * normB_global)
		{
			if (myRank == 0)
			{
				printf("Iteration: %d, Norm Reduction: %e\n",i-1, sqrt(normR_global/normB_global));// sqrt(normR_global/normB_global));
				
			}
			break;
		}

		// y =  x + D^-1 * r
		for (j = 0; j < numRow; j++)
		{
			y_local[j] = x_global[rowOrder + j] + diagInv_local[j] * r_local[j];
		}

		MPI_Allgather(y_local, numRow, MPI_DOUBLE, x_global, numRow, MPI_DOUBLE, MPI_COMM_WORLD);

		// output x1
		if (i == 1 && myRank == 0)
		{
		
			out = fopen("x1.txt", "w");
			for (j = 0; j < dim; j++)
			{
				fprintf(out,"%e\n",x_global[j]);
			}
			fclose(out);
			
		}

		// output x2
		if (i == 2 && myRank == 0)
		{
		
			out = fopen("x2.txt", "w");
			for (j = 0; j < dim; j++)
			{
				fprintf(out,"%e\n",x_global[j]);
			}
			fclose(out);
			
		}	
	}
	

	// output xfinal
	if (myRank == 0)
	{
		if (normR_global >= 10e-8 * normB_global)
		{
			printf("Jacobi did not converge: %e\n", sqrt(normR_global));
		}
		else
		{
			out = fopen("xfinal.txt", "w");
			for (j = 0; j < dim; j++)
			{
				fprintf(out,"%e\n",x_global[j]);
			}
			fclose(out);
		}
	}


	freeMatrix(A_local);
	free(x_global);
	free(y_local);
	free(diagInv_local);
	free(r_local);

	MPI_Finalize();


	return 0;
}
