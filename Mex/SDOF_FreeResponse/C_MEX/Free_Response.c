#include <stddef.h>
#include <math.h>

void Free_Response(const double* const p_w_n, const double* const p_zeta, const double* const p_x0, const double* const p_v0, const double t_vec[], double x_vec[], const size_t* const p_N)
{
	if (*p_zeta==0)
	{
		double A=sqrt(*p_x0 * *p_x0 + *p_v0 * *p_v0 / *p_w_n / *p_w_n);
		double theta=atan2(*p_x0 * *p_w_n, *p_v0);
		int n;
		for (n=0; n < (*p_N);n++)
		{
			x_vec[n]=A*sin(*p_w_n*t_vec[n]+theta);
		}
	}
	else if (*p_zeta<1)
	{
		double w_d=*p_w_n*sqrt(1- *p_zeta * *p_zeta);
		double B=sqrt(*p_x0 * *p_x0 * *p_w_n * *p_w_n+2 * *p_zeta * *p_w_n * *p_x0 * *p_v0 + *p_v0 * *p_v0) / w_d;
		double phi=atan2(*p_x0 * w_d,(*p_zeta * *p_w_n * *p_x0 + *p_v0));
		int n;
		for (n=0; n < *p_N;n++)
		{
			x_vec[n]=B * exp(-*p_zeta * *p_w_n * t_vec[n]) * sin(w_d * t_vec[n] + phi);
		}
	}
	else if (*p_zeta==1)
	{
		int n;
		for (n=0; n < *p_N;n++)
		{
			x_vec[n]=exp(- *p_w_n * t_vec[n]) * (*p_x0 + (*p_w_n * *p_x0 + *p_v0) * t_vec[n]);
		}
	}
	else if (*p_zeta > 1)
	{
		double temp=*p_w_n * sqrt(*p_zeta * *p_zeta - 1);
		double C=((temp - *p_zeta * *p_w_n) * *p_x0 - *p_v0) / 2 / temp;
		double D=((temp + *p_zeta * *p_w_n) * *p_x0 + *p_v0) / 2 / temp;
		int n;
		for (n=0; n < *p_N;n++)
		{
			x_vec[n]=exp(- *p_zeta * *p_w_n * t_vec[n]) * (C * exp(-temp * t_vec[n])+D * exp(temp * t_vec[n]));
		}
	}

	return;
}
