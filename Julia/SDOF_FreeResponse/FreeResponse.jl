function FreeResponse(w_n, zeta, x_0, v_0, t_vec)
    if zeta == 0
        A = sqrt(x_0^2 + v_0^2/w_n^2)
        theta = atan(x_0*w_n,v_0)
        x_vec=A*sin.(w_n*t_vec.+theta)
    elseif zeta<1
        w_d=w_n*sqrt(1-zeta^2)
        B=sqrt(x_0^2*w_n^2+2*zeta*w_n*x_0*v_0+v_0^2)/w_d
        phi=atan(x_0*w_d,(zeta*w_n*x_0+v_0))
        x_vec=B*exp.(-zeta*w_n*t_vec).*sin.(w_d*t_vec.+phi)
    elseif zeta==1
        x_vec=exp.(-w_n*t_vec).*(x_0.+(w_n*x_0+v_0)*t_vec)
    elseif zeta>1
        temp=w_n*sqrt(zeta^2-1)
        C=((temp-zeta*w_n)*x_0-v_0)/2/temp
        D=((temp+zeta*w_n)*x_0+v_0)/2/temp
        x_vec=exp.(-zeta*w_n*t_vec).*(C*exp.(-temp*t_vec)+D*exp.(temp*t_vec))
    end

    return x_vec
end
