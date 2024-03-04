function dxdt = odefcn_Ming_kappa2(t,x)
%% Dynamics
dxdt = zeros(2,1);
a_x=-x(1)^4-x(2)^2;
b_x=[x(1)*exp(x(2)),x(2)].';
sigma_x=sqrt(a_x^2+norm(b_x)^4);
phi=norm(b_x)^2*sqrt(1+norm(b_x)^2);
kappa=(phi-a_x*sigma_x)/(sigma_x^2+phi);
if norm(b_x)<0.01
    u=zeros(2,1);
else
    u=-(a_x+kappa*sigma_x)/(norm(b_x)^2)*b_x;
end
dxdt(1) =-x(1)^3+exp(x(2))*u(1);
dxdt(2) =-x(2)+u(2);
end