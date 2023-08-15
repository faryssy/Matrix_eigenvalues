/*

FARYSSY Mohammed
SFPN 
This program turns a randomly generated square matrix into upper hessenberg.
Args : n (matrix size)

*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>
#include <time.h>

int n; /*Matrix size */

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



void up_hessenberg(double* A)        /* Make matrix A upper Hessenberg */
{          
    double *G = malloc(n * n * sizeof(double));
    double *B = malloc(n * n * sizeof(double));

    if (G == NULL) {
        printf("Memory not allocated.\n");
    }

    if (B == NULL) {
        printf("Memory not allocated.\n");
    }
    
    for (int j = 0; j < n-2; j++) {
        for (int i = n-1; i > j+1; i--) {
            
            G = get_rotation_matrix(i-1, i, A, j);
            B = matrix_multiply(G, A);
            A = matrix_multiply(B, transpose(G));
        }
    }
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



int main(int argc, char **argv)
{
    
    n = 5; /* Default n value */

    /* Read 'n' on command line: */
        if (argc == 2)
            n = atoi(argv[1]);

    /* Initiate Matrix A */
    double *A = malloc(n * n * sizeof(double));

    srand(time(NULL));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
		    A[i * n + j] = (double)rand() / (double)RAND_MAX;
        }
	}

    /* Hessenberg*/

    //printf("Matrice A : \n");

    //display_matrix(A);

    //printf("\n \n \n");

    printf("Matrice A in hessenberg... \n");
    
    float startTime = (float)clock()/CLOCKS_PER_SEC;
    up_hessenberg(A);
    float endTime = (float)clock()/CLOCKS_PER_SEC;

    printf("Operation finished in : %f s", endTime - startTime);

    printf("\n \n \n");

    return 0;
}