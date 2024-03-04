function dxdt = odefcn_QP(t,x)
%% Dynamics -first case
dxdt = zeros(2,1);
c_x=x(1)+x(2)+2*x(1)^3*x(2);
d_x=2*x(2)^2;
Gamma_x=sqrt(c_x^2+norm(d_x)^4);
kappa=0;
if norm(d_x)<0.01
    u=0;
else
    u=-(c_x+kappa*Gamma_x)/(norm(d_x)^2)*d_x;
end
dxdt(1) =-1*x(1)-x(2);
dxdt(2) =x(1)^3+x(2)*u;
end