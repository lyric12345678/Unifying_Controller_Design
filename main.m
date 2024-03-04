clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Research Problem: Unifying Controller Design for Stabilizing Nonlinear 
%Systems with Norm-Bounded Control Inputs
%Author: Ming Li
%Date: Feb. 23. 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global m
Initial_position(:,1)=[-1,0.6];
% Initial_position(:,1)=[-5,-2.6];
t_end = 100;

%% Sontag's Universal formula Control Law
for i=1:size(Initial_position,2)
    [T,x_Sontag(:,:,i)] = ode45(@odefcn_Sontag,[0:0.1:t_end],Initial_position(:,i));
end
a_x=-x_Sontag(:,1,i).^4-x_Sontag(:,2,i).^2;
b_x=[x_Sontag(:,1,i).*exp(x_Sontag(:,2,i)),x_Sontag(:,2,i)].';
norm_b_vec=vecnorm(b_x,2,1).';
sigma_x=sqrt(a_x.^2+norm_b_vec.^4);
u_Stg_norm=abs(-(a_x+sigma_x)./(norm_b_vec.^2)).*norm_b_vec;

% figure(1)
% subplot(1,2,1)
% for i=1:1:size(Initial_position,2)
%     h_1(i)=plot(x_Sontag(:,1,i),x_Sontag(:,2,i),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('$x_1$','interpreter','latex')
% ylabel('$x_2$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,800], 'color','w')
% grid on
% subplot(1,2,2)
% u_bound=ones(size(T,1),1);
% for i=1:size(Initial_position,2)
%     h_1(i)=plot(T,abs(u_Stg_norm),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('t (s)')
% ylabel('$\|u\|$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,400], 'color','w')



%% Lin Yuandan (Universal formula with norm bounded constraints)
for i=1:size(Initial_position,2)
    [T,x_Liyuandan(:,:,i)] = ode45(@odefcn_Liyuandan,[0:0.1:t_end],Initial_position(:,i));
end
a_x=-x_Liyuandan(:,1,i).^4-x_Liyuandan(:,2,i).^2;
b_x=[x_Liyuandan(:,1,i).*exp(x_Liyuandan(:,2,i)),x_Liyuandan(:,2,i)].';
norm_b_vec=vecnorm(b_x,2,1).';
sigma_x=sqrt(a_x.^2+norm_b_vec.^4);
u_Lin_norm=abs(-(a_x+sigma_x)./(norm_b_vec.^2.*(1+sqrt(1+norm_b_vec.^2)))).*norm_b_vec;
kappa_Lin=(sigma_x-a_x.*(1+sqrt(1+norm_b_vec.^2)))./(sigma_x.*(1+sqrt(1+norm_b_vec.^2)));

% figure(2)
% subplot(1,2,1)
% for i=1:1:size(Initial_position,2)
%     h_1(i)=plot(x_Liyuandan(:,1,i),x_Liyuandan(:,2,i),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('$x_1$','interpreter','latex')
% ylabel('$x_2$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,800], 'color','w')
% grid on
% subplot(1,2,2)
% u_bound=ones(size(T,1),1);
% for i=1:size(Initial_position,2)
%     h_1(i)=plot(T,abs(u_Lin_norm),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('t (s)')
% ylabel('$\|u\|$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,400], 'color','w')



%% Ming Li(A Generalized Universal formula with norm bounded constraints) kappa_1
for i=1:size(Initial_position,2)
    [T,x_Ming_kappa1(:,:,i)] = ode45(@odefcn_Ming_kappa1,[0:0.1:t_end],Initial_position(:,i));
end
a_x=-x_Ming_kappa1(:,1,i).^4-x_Ming_kappa1(:,2,i).^2;
b_x=[x_Ming_kappa1(:,1,i).*exp(x_Ming_kappa1(:,2,i)),x_Ming_kappa1(:,2,i)].';
norm_b_vec=vecnorm(b_x,2,1).';
sigma_x=sqrt(a_x.^2+norm_b_vec.^4);
kappa_1=-2*(a_x)./sigma_x;
u_Ming_kappa1_norm=abs(-(a_x+kappa_1.*sigma_x)./(norm_b_vec.^2)).*norm_b_vec;

% figure(3)
% subplot(1,2,1)
% for i=1:1:size(Initial_position,2)
%     h_1(i)=plot(x_Ming_kappa1(:,1,i),x_Ming_kappa1(:,2,i),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('$x_1$','interpreter','latex')
% ylabel('$x_2$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,800], 'color','w')
% grid on
% subplot(1,2,2)
% u_bound=ones(size(T,1),1);
% for i=1:size(Initial_position,2)
%     h_1(i)=plot(T,u_Ming_kappa1_norm,'r-','linewidth',1.5);
%     hold on
% end
% xlabel('t (s)')
% ylabel('$\|u\|$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,400], 'color','w')




%% Ming Li (A Generalized Universal formula with norm bounded constraints) kappa_2
for i=1:size(Initial_position,2)
    [T,x_Ming_kappa2(:,:,i)] = ode45(@odefcn_Ming_kappa2,[0:0.1:t_end],Initial_position(:,i));
end
a_x=-x_Ming_kappa2(:,1,i).^4-x_Ming_kappa2(:,2,i).^2;
b_x=[x_Ming_kappa2(:,1,i).*exp(x_Ming_kappa2(:,2,i)),x_Ming_kappa2(:,2,i)].';
norm_b_vec=vecnorm(b_x,2,1).';
sigma_x=sqrt(a_x.^2+norm_b_vec.^4);
% kappa_1=(abs(a_x)-a_x)./sigma_x;
% kappa_Lin=(sigma_x-a_x.*(1+sqrt(1+norm_b_vec.^2)))./(sigma_x.*(1+sqrt(1+norm_b_vec.^2)));
kappa_2=(kappa_1+kappa_Lin)/4*3;
u_Ming_kappa2_norm=abs(-(a_x+kappa_2.*sigma_x))./norm_b_vec;

% figure(4)
% subplot(1,2,1)
% for i=1:1:size(Initial_position,2)
%     h_1(i)=plot(x_Ming_kappa2(:,1,i),x_Ming_kappa2(:,2,i),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('$x_1$','interpreter','latex')
% ylabel('$x_2$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,800], 'color','w')
% grid on
% subplot(1,2,2)
% u_bound=ones(size(T,1),1);
% for i=1:size(Initial_position,2)
%     h_1(i)=plot(T,abs(u_Ming_kappa2_norm),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('t (s)')
% ylabel('$\|u\|$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,400], 'color','w')

%% Ming Li (A Generalized Universal formula with norm bounded constraints) kappa_3
for i=1:size(Initial_position,2)
    [T,x_Ming_kappa3(:,:,i)] = ode45(@odefcn_Ming_kappa3,[0:0.1:t_end],Initial_position(:,i));
end
a_x=-x_Ming_kappa3(:,1,i).^4-x_Ming_kappa3(:,2,i).^2;
b_x=[x_Ming_kappa3(:,1,i).*exp(x_Ming_kappa3(:,2,i)),x_Ming_kappa3(:,2,i)].';
norm_b_vec=vecnorm(b_x,2,1).';
sigma_x=sqrt(a_x.^2+norm_b_vec.^4);
% kappa_1=(abs(a_x)-a_x)./sigma_x;
% kappa_Lin=(sigma_x-a_x.*(1+sqrt(1+norm_b_vec.^2)))./(sigma_x.*(1+sqrt(1+norm_b_vec.^2)));
kappa_3=(kappa_1+kappa_Lin)/3;
u_Ming_kappa3_norm=abs(-(a_x+kappa_3.*sigma_x)./(norm_b_vec).^2).*norm_b_vec;

% figure(5)
% subplot(1,2,1)
% for i=1:1:size(Initial_position,2)
%     h_1(i)=plot(x_Ming_kappa3(:,1,i),x_Ming_kappa3(:,2,i),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('$x_1$','interpreter','latex')
% ylabel('$x_2$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,800], 'color','w')
% grid on
% subplot(1,2,2)
% u_bound=ones(size(T,1),1);
% for i=1:size(Initial_position,2)
%     h_1(i)=plot(T,u_Ming_kappa3_norm,'r-','linewidth',1.5);
%     hold on
% end
% xlabel('t (s)')
% ylabel('$\|u\|$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,400], 'color','w')



%% Ming Li (A Generalized Universal formula with norm bounded constraints) optimization
for i=1:size(Initial_position,2)
    [T,x_Ming_optimization(:,:,i)] = ode45(@odefcn_Ming_optimization,[0:0.1:t_end],Initial_position(:,i));
end
a_x=-x_Ming_optimization(:,1,i).^4-x_Ming_optimization(:,2,i).^2;
b_x=[x_Ming_optimization(:,1,i).*exp(x_Ming_optimization(:,2,i)),x_Ming_optimization(:,2,i)].';
norm_b_vec=vecnorm(b_x,2,1).';
sigma_x=sqrt(a_x.^2+norm_b_vec.^4);
m=10;
for i=1:length(T)
    Verify_term(i)=m*norm(b_x(:,i))^2+sigma_x(i)^2-m*a_x(i)*norm(b_x(:,i))-m*sigma_x(i)*norm(b_x(:,i));
    if Verify_term(i)>0
        u_Ming_optimization(:,i)=-m*(a_x(i)+sigma_x(i))/(sigma_x(i)^2+m*norm(b_x(:,i))^2)*b_x(:,i);
        kappa_optimization(i)=1-(a_x(i)+sigma_x(i))*sigma_x(i)/(sigma_x(i)^2+m*norm(b_x(:,i))^2);
    end
    if Verify_term(i)<=0
        if norm(b_x(:,i))<0.01
            u_Ming_optimization(:,i)=zeros(2,1);
            kappa_optimization(i)=1;
        else
            u_Ming_optimization(:,i)=-b_x(:,i)/norm(b_x(:,i));
            kappa_optimization(i)=(norm(b_x(:,i))-a_x(i))/sigma_x(i);
        end
    end
end


% figure(6)
% subplot(1,2,1)
% for i=1:1:size(Initial_position,2)
%     h_1(i)=plot(x_Ming_optimization(:,1,i),x_Ming_optimization(:,2,i),'r-','linewidth',1.5);
%     hold on
% end
% xlabel('$x_1$','interpreter','latex')
% ylabel('$x_2$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,400], 'color','w')
% grid on
% subplot(1,2,2)
% u_bound=ones(size(T,1),1);
% for i=1:size(Initial_position,2)
%     h_1(i)=plot(T,vecnorm(u_Ming_optimization,2,1).','r-','linewidth',1.5);
%     hold on
% end
% xlabel('t (s)')
% ylabel('$\|u\|$','interpreter','latex')
% set(gca,'FontSize',23)
% set(gcf,'Position',[200,200,1000,400], 'color','w')



%% We plot out all results together
plot_step=20;
figure(1)
subplot(1,2,1)
for i=1:1:size(Initial_position,2)
    h_1(i)=plot(x_Sontag(:,1,i),x_Sontag(:,2,i),'r-','linewidth',1.5);
    hold on
    h_2(i)=plot(x_Liyuandan(:,1,i),x_Liyuandan(:,2,i),'m--','linewidth',1.5);
    h_3(i)=plot(x_Ming_kappa1(:,1,i),x_Ming_kappa1(:,2,i),'bv--','linewidth',1,'MarkerSize',4);
    h_4(i)=plot(x_Ming_kappa2(:,1,i),x_Ming_kappa2(:,2,i),'bp-','linewidth',1,'MarkerSize',4);
    h_5(i)=plot(x_Ming_kappa3(:,1,i),x_Ming_kappa3(:,2,i),'bo-','linewidth',1,'MarkerSize',4);
    h_6(i)=plot(x_Ming_optimization(:,1,i),x_Ming_optimization(:,2,i),'g+-','linewidth',1,'MarkerSize',4);
end
xlabel('$x_1$ (-)','interpreter','latex')
ylabel('$x_2$ (-)','interpreter','latex')
set(gca,'FontSize',17)
% set(gcf,'Position',[200,200,1000,400], 'color','w')
grid on
for i=1:size(Initial_position,2)
legend([h_1(i), h_2(i), h_3(i), h_4(i), h_5(i),h_6(i)], ...
    'Sontag''s Formula', 'Lin-Sontag''s Formula', 'General - $\kappa_{1}$', ...
    'General - $\kappa_{2}$', 'General - $\kappa_{3}$', 'General - Opt', ...
    'Interpreter', 'latex');
end
grid on
% axis([0,Initial_position(1),0,Initial_position(2)])
% axis([0,Initial_position(1),Initial_position(2),0])
axis([Initial_position(1),0,0,Initial_position(2)])

subplot(1,2,2)
u_opt_nom=vecnorm(u_Ming_optimization,2,1).';
u_bound=ones(size(T,1),1);
for i=1:size(Initial_position,2)
    h_1(i)=plot(T,abs(u_Stg_norm),'r-','linewidth',1.5);
    hold on
    h_2(i)=plot(T,abs(u_Lin_norm),'m--','linewidth',1.5);
    h_3(i)=plot(T,abs(u_Ming_kappa1_norm),'bv--','linewidth',1,'MarkerSize',4);
    h_4(i)=plot(T,abs(u_Ming_kappa2_norm),'bp-','linewidth',1,'MarkerSize',4);
    h_5(i)=plot(T,abs(u_Ming_kappa3_norm),'bo-','linewidth',1,'MarkerSize',4);
    h_6(i)=plot(T,u_opt_nom,'g+-','linewidth',1.5);
end
xlabel('$t$ (s)','interpreter','latex')
ylabel('$\|u\|$ (-)','interpreter','latex')
axis([0,5,0,1.5])
set(gca,'FontSize',17)
% set(gcf,'Position',[200,200,1000,800], 'color','w')
for i=1:size(Initial_position,2)
legend([h_1(i), h_2(i), h_3(i), h_4(i), h_5(i),h_6(i)], ...
    'Sontag''s Formula', 'Lin-Sontag''s Formula', 'General - $\kappa_{1}$', ...
    'General - $\kappa_{2}$', 'General - $\kappa_{3}$', 'General - Opt', ...
    'Interpreter', 'latex');
end
grid on
set(gcf,'Position',[200,200,1500,600], 'color','w')


%% Cost for each algorithm
%Lin's universal formula 
cost_Lin=u_Lin_norm.^2+m*(1-kappa_Lin).^2;
cost_kappa1=u_Ming_kappa1_norm.^2+m*(1-kappa_1).^2;
cost_kappa2=u_Ming_kappa2_norm.^2+m*(1-kappa_2).^2;
cost_kappa3=u_Ming_kappa3_norm.^2+m*(1-kappa_3).^2;
u_Ming_norm_vec=vecnorm(u_Ming_optimization,2,1).';
cost_optimization=u_Ming_norm_vec.^2+m*(1-kappa_optimization.').^2;
t_end=100;
figure(2)
for i=1:1:size(Initial_position,2)
    h_2(i)=plot(T(1:t_end),cost_Lin(1:t_end),'m--','linewidth',1.5);
    hold on
    h_3(i)=plot(T(1:t_end),cost_kappa1(1:t_end),'bv--','linewidth',1,'MarkerSize',4);
    h_4(i)=plot(T(1:t_end),cost_kappa2(1:t_end),'bp-','linewidth',1,'MarkerSize',4);
    h_5(i)=plot(T(1:t_end),cost_kappa3(1:t_end),'bo-','linewidth',1,'MarkerSize',4);
    h_6(i)=plot(T(1:t_end),cost_optimization(1:t_end),'g+-','linewidth',1,'MarkerSize',4);
end
xlabel('$t$ (s)','interpreter','latex')
ylabel('Cost (-)')
set(gca,'FontSize',17)
% set(gcf,'Position',[200,200,1000,800], 'color','w')
for i=1:size(Initial_position,2)
legend([h_2(i), h_3(i), h_4(i), h_5(i),h_6(i)], ...
    'Lin-Sontag''s Formula', 'General - $\kappa_{1}$', ...
    'General - $\kappa_{2}$', 'General - $\kappa_{3}$', 'General - Opt', ...
    'Interpreter', 'latex');
end
grid on
set(gcf,'Position',[200,200,750,600], 'color','w')