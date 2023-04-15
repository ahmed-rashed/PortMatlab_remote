using Plots
include("FreeResponse.jl")

x_0=-1
v_0=0
w_n=1
f_n=w_n/2pi
T_n=1/f_n

zeta_vec=[0,.1,.2,.4,1/sqrt(2),1,2]
legend_str = ["\$\\zeta = $i \$" for i in union(0:0.1:0.2,[0.4,1/sqrt(2),1,2])]
legend_str[5] = "\$\\zeta = 1/\\sqrt{2}\$"
t_vec= 0:2*T_n/(500-1):2*T_n

plot()
for (n, zeta) in enumerate(zeta_vec)
    x_vec=FreeResponse(w_n,zeta,x_0,v_0,t_vec)
    plot!(t_vec/T_n,x_vec,label = legend_str[n])
end
p = plot!(legend=:bottomright,grid = true)

title!("\$x(t)\$ for \$\\omega_{n}=1\$, \$x_{0}=-1\$ and \$v_{0}=0\$")
xlabel!("\$t/T_{n}\\qquad,:T_{n}=1/f_{n}=2\\pi/\\omega_{n}\$")

display(p)

png(p,"SDOF_Julia")
Plots.pdf(p, "SDOF_Julia")
Plots.svg(p, "SDOF_Julia")