#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double normHill_mod(double x,double n,double k)
{
	if(x<0) x=0.0;
	else if(x>1) x=1.0;
	if(x==0.0 && k==0.0) return(1.0);
	return(pow(x,exp(2*k))/(pow(x,exp(2*k))+pow(k,exp(2*k)))*(1+pow(k,exp(2*k)))+0*n);
}
