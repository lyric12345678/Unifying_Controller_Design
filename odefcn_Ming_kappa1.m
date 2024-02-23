function dxdt = odefcn_Ming_kappa1(t,x)
%% Dynamics
dxdt = zeros(2,1);
a_x=-x(1)^4-x(2)^2;
b_x=[x(1)*exp(x(2)),x(2)].';
sigma_x=sqrt(a_x^2+norm(b_x)^4);
kappa=(0.1*norm(a_x)-a_x)/sigma_x;
% if (-a_x)/sigma_x>0
%     kappa=(-a_x)/sigma_x;
% else
%     kappa=0;
% end
% kappa=1/10*log(exp((-10*a_x)/sigma_x)+exp(-10*0));
if norm(b_x)<0.01
    u=zeros(2,1);
else
    u=-(a_x+kappa*sigma_x)/(norm(b_x)^2)*b_x;
end
dxdt(1) =-x(1)^3+exp(x(2))*u(1);
dxdt(2) =-x(2)+u(2);
end