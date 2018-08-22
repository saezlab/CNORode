#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double normHill(double x,double n,double k)
{
	if(x<0) x=0.0;
	else if(x>1) x=1.0;
	if(x==0.0 && k==0.0) return(1.0);
	return(pow(x,n)/(pow(x,n)+pow(k,n))*(1+pow(k,n)));
}
