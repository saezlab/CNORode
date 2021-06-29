#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double FG_transfer_function(double x,double n,double k)
{
	return(1-pow(1-x,n)/(pow(1-x,n)+pow(k,n))*(1+pow(k,n)));
}
