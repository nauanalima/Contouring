#include <stdio.h>
#include <math.h>

#define n 2
#define np 100
#define h 0.5
#define x0 0
#define x1 1

typedef double (*systemfunc)();

double f0 (double x, double *y){
	return (y[1]);
}

double f2 (double x, double *y) {
	return (4*y[0]);
}

double f1 (double x, double *y) {
	return (f2(x,y)-4*x);
}

double shooting (systemfunc func[], double x, double *y) {
	double k1[n], k2[n], k3[n], k4[n], temp[n];
	int i;
	
	for(i=0; i<n; i++) {
		k1[i] = func[i](x,y);
	}
	for(i=0; i<n; i++) {
		temp[i] = y[i]+h*k1[i]/2.;
	}
	for(i=0; i<n; i++) {
		k2[i] = func[i](x+h/2.,temp);
	}
	for(i=0; i<n; i++) {
		temp[i] = y[i]+h*k2[i]/2.;
	}
	for(i=0; i<n; i++) {
		k3[i] = func[i](x+h/2.,temp);
	}
	for(i=0; i<n; i++) {
		temp[i] = y[i]+h*k3[i];
	}
	for(i=0; i<n; i++) {
		k4[i] = func[i](x+h,temp);	
	}
	for(i=0; i<n; i++) {
		y[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6.;
	}
}

int main() {
	int i;
	double x, y, y1b, y2b, b = 2;
	double y1[n]={0,0}, y2[n]={0,1};
	double yp1[np+1], yp2[np+1];
	systemfunc equations1[n]={f0,f1};
	systemfunc equations2[n]={f0,f2};

	FILE *f;
	f = fopen("prob1a.dat", "w+");

	i=0;
	yp1[0] = y1[0];
	yp2[0] = y2[0];
	for(x=x0+h; x<=x1; x+=h) {	
		yp1[i]=shooting(equations1,x,y1);
		yp2[i]=shooting(equations2,x,y2);
		i++;
	}

	i=0;
	y1b = yp1[i];
	y2b = yp2[i];
	for( x=x0; x<=x1; x+=h) {	
		y = yp1[i]+((b-y1b)/y2b)*yp2[i]; 
		fprintf(f, "%lf\t%lf\n", x, y);
		printf("%lf\t%lf\n", x, y);
		i++;
	}
	fclose(f);
	return 0;
}
