function dxdt = odefcn_Ming_optimization(t,x)
global m
%% Dynamics
dxdt = zeros(2,1);
a_x=-x(1)^4-x(2)^2;
b_x=[x(1)*exp(x(2)),x(2)].';
sigma_x=sqrt(a_x^2+norm(b_x)^4);
m_bound=a_x*sigma_x/(norm(b_x)^2);
m=1000;
Verify_term=m*norm(b_x)^2+sigma_x^2-m*a_x*norm(b_x)-m*sigma_x*norm(b_x);
if Verify_term>0
    u=-m*(a_x+sigma_x)/(sigma_x^2+m*norm(b_x)^2)*b_x;
end
if Verify_term<=0
    if norm(b_x)<0.01
        u=zeros(2,1);
    else
        u=-b_x/abs(b_x);
    end
end
dxdt(1) =-x(1)^3+exp(x(2))*u(1);
dxdt(2) =-x(2)+u(2);
end