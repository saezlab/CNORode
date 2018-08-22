#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double FG_transfer_function(double x,double n,double k)
{
	if(x<0) x=0.0;
	else if(x>1) x=1.0;
	if(x==1.0 && k==0.0) return(0.0);
	return(1-pow(1-x,n)/(pow(1-x,n)+pow(k,n))*(1+pow(k,n)));
}
