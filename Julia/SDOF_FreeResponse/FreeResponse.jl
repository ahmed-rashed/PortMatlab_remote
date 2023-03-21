function FreeResponse(ωₙ, ζ, x₀, v₀, t_vec)
    if ζ == 0
        A = √(x₀^2 + v₀^2/ωₙ^2)
        θ = atan(x₀*ωₙ,v₀)
        x_vec=A*sin.(ωₙ*t_vec.+θ)
    elseif ζ<1
        w_d=ωₙ*sqrt(1-ζ^2)
        B=sqrt(x₀^2*ωₙ^2+2*ζ*ωₙ*x₀*v₀+v₀^2)/w_d
        ϕ=atan(x₀*w_d,(ζ*ωₙ*x₀+v₀))
        x_vec=B*exp.(-ζ*ωₙ*t_vec).*sin.(w_d*t_vec.+ϕ)
    elseif ζ==1
        x_vec=exp.(-ωₙ*t_vec).*(x₀.+(ωₙ*x₀+v₀)*t_vec)
    elseif ζ>1
        temp=ωₙ*sqrt(ζ^2-1)
        C=((temp-ζ*ωₙ)*x₀-v₀)/2/temp
        D=((temp+ζ*ωₙ)*x₀+v₀)/2/temp
        x_vec=exp.(-ζ*ωₙ*t_vec).*(C*exp.(-temp*t_vec)+D*exp.(temp*t_vec))
    end

    return x_vec
end
