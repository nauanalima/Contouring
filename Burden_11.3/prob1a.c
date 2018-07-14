#include <stdio.h>
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

void exchangelines (double *line1, double *line2, int row) {
	double temp;
	int i;

	for (i=0; i<=row; i++) {
		temp = line1[i];
		line1[i] = line2[i];
		line2[i] = temp;
	}
}

double uppertriangular (double **matrix, int row, int col) {
	int i, j, k, k1, test, pivot;
	double m, temp;

	for (k=0; k<row; k++) {
		for (j=k; j<row; j++) {
			test = -1;
			if (fabs(matrix[j][k]) > matrix[k][k]) {
				test = j;
			}
			if (test!=-1) {
				exchangelines(matrix[k], matrix[test], row);
			}
		}
		for (i=k+1; i<row; i++) {
			m = matrix[i][k]/matrix[k][k];
			for (j=k; j<=col; j++)
				matrix[i][j] = matrix[i][j] - m*matrix[k][j];
		}
	}
	return **matrix;		
}

double *reversesub (double **matrix, const int dim) {
	int j, k, n;
	double sum, *roots = (double*)malloc(dim*sizeof(double));
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
	int i, j;
    double h = 0.5;
}
