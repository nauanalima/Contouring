#include <stdio.h>
#include <math.h>

#define n 2
#define np 100

typedef double (*systemfunc)();

double f0 (double x, double *y){
	return (y[1]);
}

double f2 (double x, double *y) {
	return (100*y[0]);
}

double f1 (double x, double *y) {
	return (f2(x,y));
}

double shooting (systemfunc func[], double x, double *y, double h) {
	double k1[n], k2[n], k3[n], k4[n], temp[n];
	int i;
	
	for(i=0; i<n; i++) {
		k1[i] = func[i](x,y);
		temp[i] = y[i]+h*k1[i]/2.;
	}
	for(i=0; i<n; i++) {
		k2[i] = func[i](x+h/2.,temp);
		temp[i] = y[i]+h*k2[i]/2.;
	}
	for(i=0; i<n; i++) {
		k3[i] = func[i](x+h/2.,temp);
		temp[i] = y[i]+h*k3[i];
	}
	for(i=0; i<n; i++) {
		k4[i] = func[i](x+h,temp);	
		y[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6.;
	}
}

int main() {
	int i;
	double x, h, min=0, max=1, y, y1b, y2b, b=exp(1e-10);
	double y1[n]={1,0}, y2[n]={0,1};
	double yp1[np+1], yp2[np+1];
	systemfunc equations1[n]={f0,f1};
	systemfunc equations2[n]={f0,f2};
	
	h=0.05;

	i=0;
	yp1[0] = y1[0];
	yp2[0] = y2[0];
	for(x=min+h; x<=max; x+=h) {	
		yp1[i]=shooting(equations1,x,y1,h);
		yp2[i]=shooting(equations2,x,y2,h);
		i++;
	}

	i=0;
	y1b = yp1[i];
	y2b = yp2[i];
	for( x=min; x<=max; x+=h) {	
		y = yp1[i]+((b-y1b)/y2b)*yp2[i]; 
		printf("%lf\t%lf\n", x, y);
		i++;
	}
	
	return 0;
}
