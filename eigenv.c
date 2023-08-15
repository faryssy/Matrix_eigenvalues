/*

FARYSSY Mohammed
SFPN 
This program computes approximatively the eigenvalues of a random matrix A using givens' algorithm.
The accuracy of the results is verified using the trace property, i.e

trace = sum of diagonal coeff = sum of eigenvalues;
Args : n (matrix size)
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>
#include <time.h>

int n; /*Matrix size */

double trace(double*A)
{
    double t = 0;
    for (int i = 0; i < n; i++)
    {
     t += A[i * n + i];   
    }
    return t;
}

double* matrix_multiply(double *A, double *B)         /* Matrix multiplication */
{

    double *C = malloc(n * n * sizeof(double));
    if (C == NULL) {
        printf("Memory not allocated.\n");
    }
    
    for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
                C[i * n + j] = 0;
				for (int k = 0; k < n; k++)
					C[i * n + j] += A[i * n + k] * B[k * n + j]; 
            }

    return C;
}

double* transpose(double* A)               /* Transpose matrix A */
{
    double *tA = malloc(n * n * sizeof(double));
    if (tA == NULL) {
        printf("Memory not allocated.\n");
    }

    for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
                tA[i * n +j] = A[j * n + i];

    return tA;
}


void display_matrix(double* A)              /* Display the matrix A */
{
     for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
		    printf("%.2f ", A[i * n + j]);
        }

        printf("\n");
	}
}

double* get_rotation_matrix(int i1, int i2, double* A, int col)  /* Generate rotation matrix G_{i1,i2} based on matrix A*/
{
    double *G = malloc(n * n * sizeof(double));
    if (G == NULL) {
        printf("Memory not allocated.\n");
    }

    double denum = sqrt(pow(A[(i1)*n + col],2) + pow(A[(i2)*n + col], 2));
    double c = A[(i1)* n + col]/denum;
    double s = A[(i2)* n + col]/denum;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {

            if (j == i) {
                if (i == i1 || i == i2)
                    G[i * n + j] = c;

                else
                    G[i * n + j] = 1;
            }

            else if (i == i1 && j == i2) {
                G[i * n + j] = s;

            } else if (i == i2 && j == i1) {
                G[i * n + j] = -s;
                
            } else { 
                G[i * n + j] = 0;
            }    
        }
    }

    return G;
}

double* givens(double* A)        /* Return Q matrix in the QR decomposition using givens */
{          
    double *G = malloc(n * n * sizeof(double));
    double *Q = malloc(n * n * sizeof(double));


    if (G == NULL) {
        printf("Memory not allocated.\n");
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {   
            if (i == j)
                Q[i * n + j] = 1;
            else
                Q[i * n + j] = 0;
        }
    }
    
    for (int j = 0; j < n; j++) 
    {
        for (int i = j+1; i < n; i++) 
        {
            G = get_rotation_matrix(j, i, A, j);
            Q = matrix_multiply(Q, transpose(G));
        }
    }

    return Q;
}

double subdiag_mean(double* A) 
{
    double subdiag_mean = 0;

    for (int i = 1; i < n; i++)
    {
       subdiag_mean += fabs(A[i * n + (i-1)]); 
    }
    
    
    return subdiag_mean/(n-1);
}

double* eigenv(double* A, double* E, double threshold)        /* Compute eigenvalues using givens' algorithm */
{
           
    double *Q = malloc(n * n * sizeof(double));

    if (Q == NULL)
    {
        printf("Memory not allocated.\n");
    }

    Q = givens(A);
    double s = subdiag_mean(A);

    if (s < threshold) {                 // If the coefficients on the first subdiagonal are under the threshold, retrieve the diagonal coefficients as eigenvalues.
        for (int i = 0; i < n; i++){
            E[i] = A[i * n + i];
        }
        return E;
   }
    
    A = matrix_multiply(transpose(Q), A);
    A = matrix_multiply(A, Q);
    
    return eigenv(A, E, threshold);
   
}



void display_eig(double * E)
{
    for (int i = 0; i < n; i++)
    {
        printf("%.2f  ", E[i]);
    }

    printf("\n");
}



int main(int argc, char **argv)
{
    
    n = 5; /* Default n value */

    /* Read 'n' on command line: */
        if (argc == 2)
            n = atoi(argv[1]);

    /* Initiate Matrix A and eigenvalues list E */
    double *A = malloc(n * n * sizeof(double));
    double *E = malloc(n * sizeof(double));
    
    srand(time(NULL));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {   
		    A[i * n + j] = (double)rand() / (double)RAND_MAX;
        }

        E[i] = 0;
	}
 
    double sl;
    

    printf("\n Matrice A : \n");

    display_matrix(A);

    printf("\n Computing eigenvalues... \n");
    
    float startTime = (float)clock()/CLOCKS_PER_SEC;
    E = eigenv(A, E, 0.1);
    printf("The result is : \n \n");
    float endTime = (float)clock()/CLOCKS_PER_SEC;

    for (int i = 0; i < n; i++)
    {
        sl+= E[i];
    }

    display_eig(E);

    printf("\nThis is : %f far from the accurate values.", fabs(sl/n - trace(A)/n));
    printf("\n");
    
    printf("\n Operation finished in : %f s", endTime - startTime);

    return 0;
}