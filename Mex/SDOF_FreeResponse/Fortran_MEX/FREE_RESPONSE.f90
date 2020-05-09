FUNCTION FREE_RESPONSE(w_n,zeta,x0,v0,t_vec,N)
IMPLICIT NONE
INTEGER(4), INTENT(IN) :: N
REAL(8) w_n,zeta,x0,v0,t_vec(N),FREE_RESPONSE(N)
REAL(8) A,theta,w_d,B,phi,temp,C,D


IF (zeta==0) THEN
    A=DSQRT(x0**2+v0**2/w_n**2)
    theta=DATAN2(x0*w_n,v0)
    FREE_RESPONSE=A*DSIN(w_n*t_vec+theta)
ELSE IF (zeta<1) THEN
    w_d=w_n*DSQRT(1-zeta**2)
    B=DSQRT(x0**2*w_n**2+2*zeta*w_n*x0*v0+v0**2)/w_d
    phi=DATAN2(x0*w_d,(zeta*w_n*x0+v0))
    FREE_RESPONSE=B*DEXP(-zeta*w_n*t_vec)*DSIN(w_d*t_vec+phi)
ELSE IF (zeta==1) THEN
    FREE_RESPONSE=DEXP(-w_n*t_vec)*(x0+(w_n*x0+v0)*t_vec)
ELSE IF (zeta>1)  THEN
    temp=w_n*DSQRT(zeta**2-1)
    C=((temp-zeta*w_n)*x0-v0)/2/temp
    D=((temp+zeta*w_n)*x0+v0)/2/temp
    FREE_RESPONSE=DEXP(-zeta*w_n*t_vec)*(C*DEXP(-temp*t_vec)+D*DEXP(temp*t_vec))
END IF

END FUNCTION FREE_RESPONSE
