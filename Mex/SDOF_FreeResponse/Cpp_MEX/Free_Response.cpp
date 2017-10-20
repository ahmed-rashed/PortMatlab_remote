#include <valarray>
#include <math.h>

std::valarray<double> Free_Response(const double& w_n, const double& zeta, const double& x0, const double& v0, const std::valarray<double>& t_vec)
{
	if (zeta==0)
	{
		double A=sqrt(x0*x0+v0*v0/w_n/w_n);
		double theta=atan2(x0*w_n,v0);
		return A*sin(w_n*t_vec+theta);
	}
	else if (zeta<1)
	{
		double w_d=w_n*sqrt(1-zeta*zeta);
		double B=sqrt(x0*x0*w_n*w_n+2*zeta*w_n*x0*v0+v0*v0)/w_d;
		double phi=atan2(x0*w_d,(zeta*w_n*x0+v0));
		return B*exp(-zeta*w_n*t_vec)*sin(w_d*t_vec+phi);
	}
	else if (zeta==1)
	{
		return exp(-w_n*t_vec)*(x0+(w_n*x0+v0)*t_vec);
	}
	else if (zeta>1)
	{
		double temp=w_n*sqrt(zeta*zeta-1);
		double C=((temp-zeta*w_n)*x0-v0)/2/temp;
		double D=((temp+zeta*w_n)*x0+v0)/2/temp;
		return exp(-zeta*w_n*t_vec)*(C*exp(-temp*t_vec)+D*exp(temp*t_vec));
	}

	return std::valarray<double>();
}