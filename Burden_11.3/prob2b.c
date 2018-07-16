#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 1
#define pi 3.14159265359
#define x0 0
#define x1 pi/2
#define y0 -0.3
#define y1 -0.1

double p (double x) {
	return (0);
}

double q (double x) {
	return (2);
}

double r (double x) {
	return (cos(x));
}

double y (double x) {
    return (-(sin(x)+3*cos(x)/10));
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

double **creatematrix (double **matrix, double x[], int h) {
    int i, j;
	
	for(i=0;i<n;i++) {
		for(j=0;j<(n+1);j++)
			matrix[i][j] = 0;
    }

	for(i=1;i<n-1;i++) {
		matrix[i][n] = -1*pow(h,2)*r(x[i+1]);
		matrix[i][i] = 2 + pow(h,2)*q(x[i+1]);
		matrix[i][i-1] = -1 - (h/2.0)*p(x[i+1]);
		matrix[i][i+1] = -1 + (h/2.0)*p(x[i+1]);		
	}
	matrix[0][0] = 2 + pow(h,2)*q(x[1]);
	matrix[0][1] = -1 + (h/2.0)*p(x[1]);
    matrix[0][n] = -1*pow(h,2)*r(x[2]) + (1+(h/2)*p(x[1]))*y0;
	matrix[n-1][n-1] = 2 + pow(h,2)*q(x[1]);
	matrix[n-1][n-2] = -1 -(h/2.0)*p(x[n]);
	matrix[n-1][n] = -1*pow(h,2)*r(x[n]) + (1 - (h/2.0)*p(x[n]))*y1;
	
    return matrix;
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
	int i, j, row = n, col = n+1, steps;
    double h = pi/8, x[n+2], X = x0;
    double **matrix, *roots;

	FILE *f;
	f = fopen("prob2b.dat", "w+");

    matrix = malloc(row* sizeof(double*));
	for( i = 0; i < col; i++ )
        matrix[i] = malloc(col*sizeof(double));

    x[0] = x0;

    for(i=1; i<=n; i++)	
		x[i] = x0 + i*h;
    x[n+1] = x1;

    creatematrix(matrix,x,h);
    print(matrix,row,col);
    steps = uppertriangular(matrix,row,col);
    print(matrix,row,col);
    roots = reversesub(matrix,row);

	printf("%lf\t%lf\n", X, roots[i]);
	fprintf(f,"%lf\t%lf\n",X, y0);
    for(i=0;i<n;i++) {
        printf("%lf\t%lf\n", X+h, roots[i]);
		fprintf(f,"%lf\t%lf\n",X+h, roots[i]);
		X+=h;
	}
	printf("%lf\t%lf\n", X+h, y1);
	fprintf(f,"%lf\t%lf\n",X+h, y1);

	printf("\n\n\n\n\n");
}
