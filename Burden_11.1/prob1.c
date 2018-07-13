#include <stdio.h>

double f(double x, double y){
	return (4*(y-x));
}

double shooting (double x, double y, double h) {
	double k1, k2, k3, k4;
	
	k1=f(x,y);
	k2=f(x+0.5*h,y+0.5*h*k1);
	k3=f(x+0.5*h,y+0.5*h*k2);
	k4=f(x+h,y+h*k3);
	
	y=y+(h/6)*(k1+2*k2+2*k3+k4);
	
	return y;
}

int main(int argc, char **argv) {
	int N=2;
	double x, h, y;
	double a=0, b=1;
	
	h=(b-a)/N;
	y=1;
	
	return 0;
}