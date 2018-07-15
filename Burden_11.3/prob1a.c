#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 3

double p (double x) {
	return (0);
}

double q (double x) {
	return (4);
}

double r (double x) {
	return (-4*x);
}

double y (double x) {
    return (exp(2)*(exp(2*x)-exp(-2*x))/(exp(4)-1)+x);
}

void print(double **matrix, int row, int col){
	int i, j;
	printf("\nMatriz: \n");
	for(i=0;i<row;i++) {
		for(j=0;j<col;j++) { 
			printf("%5.2lf\t",matrix[i][j]);
		}	  
		printf("\n");
	}
	printf("\n");
}
double **creatematrix (double **matrix, int row, int col, int *x, int h) {
    int i;

    for(i=0; i<row; i++) {	
		if(i==0) {
			matrix[i][i] = 2 + pow(h,2)* q(x[i+1]);
			matrix[i][i+1] = - 1 + (h/2.)* p(x[i+1]);
		}
		else {
            if(i==n-1) {
            	matrix[i][i-1] = - 1 - (h/2.)* p(x[i+1]);
                matrix[i][i] = 2 + pow(h, 2)* q(x[i+1]);
            }
            else {
            	matrix[i][i-1] = - 1 - (h/2.)* p(x[i+1]);
                matrix[i][i] = 2 + pow(h, 2)* q(x[i+1]);
                matrix[i][i+1] = - 1 + (h/2.)* p(x[i+1]);
            }
		}
	}
	matrix[0][n] = -pow(h,2)*r(x[1]) + (1+(h/2.)*p(x[1]))*y(x[0]);
    matrix[n-1][n] = -pow(h,2)*r(x[n]) + (1-(h/2.)*p(x[n]))*y(x[n+1]);
}

double jacobi(double **matrix, int row, int col, double *x1, double *x2) {

	int i,j;
	double sum, s1, s2, error;

	for(i=0; i<row; i++) {	
		sum = 0;
		
		for(j=0; j<col-1; j++) {	
			if( j!= i )
				sum += -matrix[i][j]*x1[j];
		}

		x2[i] = (sum + matrix[i][col-1])/ matrix[i][i];
	}

	for(i=0; i<row; i++) {	

		s1 += pow(x2[i]-x1[i], 2);
		s2 += pow(x2[i], 2);
	
	}
	
	error = sqrt(s1)/sqrt(s2);
	
	return error;
}

void exchangelines (double *line1, double *line2, int row) {
	double temp;
	int i;

	for (i=0; i<=row; i++) {
		temp = line1[i];
		line1[i] = line2[i];
		line2[i] = temp;
	}
}

int uppertriangular (double **matrix, int row, int col) {
	int i, j, k, k1, test, pivot;
	double m, temp;

	int steps = 0; 

	for (k=0; k<row; k++) {
		for (j=k; j<row; j++) {
			test = -1;
			if (fabs(matrix[j][k]) > matrix[k][k]) {
				test = j;
			}
			if (test!=-1) {
				exchangelines(matrix[k], matrix[test], row);
				steps++;
			}
		}
		for (i=k+1; i<row; i++) {
			m = matrix[i][k]/matrix[k][k];
			for (j=k; j<=col; j++)
				matrix[i][j] = matrix[i][j] - m*matrix[k][j];
		}
	}
	return steps;		
}

double *reversesub (double **matrix, const int dim) {
	int j, k;
	double sum;
    double *roots = malloc(dim*sizeof(double));
	roots[dim-1] = matrix[dim-1][dim]/matrix[dim-1][dim-1];
	k = dim-2;
	while (k>=0) {
		sum = matrix[k][dim];
		for (j=k+1; j<dim; j++) 
			sum = sum - matrix[k][j]*roots[j];
		roots[k] = sum/matrix[k][k];
		k--;
	}
	return roots;
}

int main() {
	int i, j, row = n, col = n+1;
    double h = 0.5, x[n+2], w[n+2], x0, x1, y0, y1;
    double **matrix;

    matrix = malloc(row* sizeof(double*));
	for( i = 0; i < col; i++ )
        matrix[i] = malloc(col*sizeof(double));

    x[0] = x0;

    for(i=1; i<=n; i++)	
		x[i] = x0 + i*h;
    x[n+1] = x1;

}
