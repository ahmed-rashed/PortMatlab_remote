clc
clearvars
close all

x0=-1;
v0=0;
w_n=1;
f_n=w_n/2/pi;
T_n=1/f_n;

zeta_vec=[0,.1,.2,.4,1/sqrt(2),1,2];
legend_string={'$\zeta = 0$','$\zeta = 0.1$','$\zeta = 0.2$','$\zeta = 0.4$','$\zeta=1/\sqrt{2}$','$\zeta = 1$','$\zeta = 2$'};

t_vec=linspace(0,2*T_n,500);

figure
hold all
for n=1:length(zeta_vec)
    x_vec=Free_Response(w_n,zeta_vec(n),x0,v0,t_vec);
    x_vec1=C_MEX(w_n,zeta_vec(n),x0,v0,t_vec);clear C_MEX;max(abs(x_vec1-x_vec)),x_vec=x_vec1;
    x_vec2=Cpp_MEX(w_n,zeta_vec(n),x0,v0,t_vec);clear Cpp_MEX;max(abs(x_vec2-x_vec)),x_vec=x_vec2;
    x_vec3=Fortran_MEX(w_n,zeta_vec(n),x0,v0,t_vec);clear Fortran_MEX;max(abs(x_vec3-x_vec)),x_vec=x_vec3;
    plot(t_vec/T_n,x_vec)
end

title('$x(t)$ for $\omega_{n}=1$, $x_{0}=-1$ and $v_{0}=0$','interpreter','latex');
xlabel('$t/T_{n}\qquad,:T_{n}=1/f_{n}=2\pi/\omega_{n}$','interpreter','latex');

legend(legend_string,'interpreter','latex','Location','SouthEast');
grid on

export_figure(gcf,'',{'SDOF_1'})