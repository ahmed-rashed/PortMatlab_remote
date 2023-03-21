using Plots
include("FreeResponse.jl")

x₀=-1
v₀=0
ωₙ=1
fₙ=ωₙ/2/π
Tₙ=1/fₙ

zeta_vec=[0,.1,.2,.4,1/√2,1,2]
legend_str = ["\$\\zeta = $i \$" for i in ∪(0:0.1:0.2,[0.4,1/√2,1,2])]
legend_str[5] = "\$\\zeta = 1/\\sqrt{2}\$"
t_vec= 0:2*Tₙ/(500-1):2*Tₙ

plot()
for (n, ζ) in enumerate(zeta_vec)
    x_vec=FreeResponse(ωₙ,ζ,x₀,v₀,t_vec)
    plot!(t_vec/Tₙ,x_vec,label = legend_str[n])
end
p = plot!(legend=:bottomright,grid = true)

title!("\$x(t)\$ for \$\\omega_{n}=1\$, \$x_{0}=-1\$ and \$v_{0}=0\$")
xlabel!("\$t/T_{n}\\qquad,:T_{n}=1/f_{n}=2\\pi/\\omega_{n}\$")

display(p)

png(p,"SDOF_Julia")
Plots.pdf(p, "SDOF_Julia")
Plots.svg(p, "SDOF_Julia")