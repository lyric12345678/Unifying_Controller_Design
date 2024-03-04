function dxdt = odefcn_Sontag(t,x)
%% Dynamics -first case
dxdt = zeros(2,1);
a_x=-x(1)^4-x(2)^2;
b_x=[x(1)*exp(x(2)),x(2)].';
sigma_x=sqrt(a_x^2+norm(b_x)^4);
if norm(b_x)<0.01
    u=zeros(2,1);
else
    u=-(a_x+sigma_x)/(norm(b_x)^2)*b_x;
end
dxdt(1) =-x(1)^3+exp(x(2))*u(1);
dxdt(2) =-x(2)+u(2);


% %% Dynamics -first case
% dxdt = zeros(2,1);
% a_x=0;
% b_x=[x(1),x(2)].';
% sigma_x=sqrt(a_x^2+norm(b_x)^4);
% if norm(b_x)<0.01
%     u=zeros(2,1);
% else
%     u=-(a_x+sigma_x)/(norm(b_x)^2)*b_x;
% end
% dxdt(1) =x(2)+u(1);
% dxdt(2) =-x(1)+u(2);


% %% Dynamics -second case
% dxdt = zeros(2,1);
% V_x=1/2*(x(1)^2+x(2)^2);
% a_x=x(1)*x(2)-x(2)*x(1)^3+2*V_x;
% b_x=[x(1)*sin(x(1)),x(2)].';
% sigma_x=sqrt(a_x^2+norm(b_x)^4);
% if norm(b_x)<0.01
%     u=zeros(2,1);
% else
%     u=-(a_x+sigma_x)/(norm(b_x)^2)*b_x;
% end
% dxdt(1) =x(2)+sin(x(1))*u(1);
% dxdt(2) =-x(1)^3+u(2);
end