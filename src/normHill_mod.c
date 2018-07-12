#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double normHill_mod(double x,double n,double k)
{
	return(pow(x,exp(2*k))/(pow(x,exp(2*k))+pow(k,exp(2*k)))*(1+pow(k,exp(2*k)))+0*n);
}
