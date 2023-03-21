using Plots
function Free_Response(ωₙ, ζ, x₀, v₀, t_vec)
    if ζ == 0
        A = √(x₀^2 + v₀^2/ωₙ^2);
        θ = atan(x₀*ωₙ,v₀);
        x_vec=A*sin(ωₙ*t_vec+θ);
    elseif ζ<1
        w_d=ωₙ*sqrt(1-ζ^2);
        B=sqrt(x₀^2*ωₙ^2+2*ζ*ωₙ*x₀*v₀+v₀^2)/w_d;
        ϕ=atan(x₀*w_d,(ζ*ωₙ*x₀+v₀));
        x_vec=B*exp(-ζ*ωₙ*t_vec).*sin(w_d*t_vec+ϕ);
    elseif ζ==1
        x_vec=exp(-ωₙ*t_vec).*(x₀+(ωₙ*x₀+v₀)*t_vec);
    elseif ζ>1
        temp=ωₙ*sqrt(ζ^2-1);
        C=((temp-ζ*ωₙ)*x₀-v₀)/2/temp;
        D=((temp+ζ*ωₙ)*x₀+v₀)/2/temp;
        x_vec=exp(-ζ*ωₙ*t_vec).*(C*exp(-temp*t_vec)+D*exp(temp*t_vec));
    end
end

x₀=-1;
v₀=0;
ωₙ=1;
fₙ=ωₙ/2/π;
Tₙ=1/fₙ;

zeta_vec=[0,.1,.2,.4,1/√2,1,2];
legend_str = ["\$\\zeta = $i \$" for i in ∪(0:0.1:0.2,[0.4,1/√2,1,2])]
legend_str[5] = "\$\\zeta = 1/\\sqrt{2}\$"
t_vec= 0:2*Tₙ/(500-1):2*Tₙ

plot()
for (n, ζ) in enumerate(zeta_vec)
    x_vec=Free_Response.(ωₙ,ζ,x₀,v₀,t_vec);
    plot!(t_vec/Tₙ,x_vec,label = legend_str[n])
end
p = plot!(legend=:bottomright,grid = true)

title!("\$x(t)\$ for \$\\omega_{n}=1\$, \$x_{0}=-1\$ and \$v_{0}=0\$")
xlabel!("\$t/T_{n}\\qquad,:T_{n}=1/f_{n}=2\\pi/\\omega_{n}\$")

png(p,"SDOF_1")
Plots.pdf(p, "SDOF_1")
Plots.svg(p, "SDOF_1")