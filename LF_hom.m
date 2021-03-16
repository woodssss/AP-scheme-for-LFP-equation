function LF_hom(N,alpha,epsi,L,l_lim,IC)
%% fixed 
close all
warning off;
ds = pi/N;
s=(pi/N/2:pi/N:2*pi-pi/2/N);
switch IC
    case 1
        f_0 = exp(-(L*cot(s)).^2);
        disp('exponential decay IC')
    case 2
        f_0=1./(1+(L*cot(s)).^2).^((3+alpha)/2);
        disp('Power law IC with -(3+alpha)')
    case 3
        f_0=1./(1+(L*cot(s)).^2).^((2+alpha)/2);
        disp('Power law IC with -(2+alpha)')
    case 4
        f_0=1./(1+(L*cot(s)).^2).^((1+alpha)/2);
        disp('Power law IC with -(1+alpha)')
    case 5
        f_0=1./(1+(L*cot(s)).^2).^((1+alpha/2)/2);
        disp('Power law IC with -(1+alpha)')
    otherwise
        error('Case not available, choose from {1,2,3,4,5}')
end

%fl_exact = 2^alpha*gamma(0.5+alpha/2)*gamma(0.5-alpha/2)^(-1)*(1+x.^2).^(-(1+alpha)/2);
T=40;
dt=0.01;
%dt=0.01;
t=0;
f=f_0';
xs=L*cot(s);
%% Implicit evolution matrix
[P]=LFP_homo_mat(N,alpha,L,l_lim,dt,epsi);
weight=L./((sin(s(1:N))).^2);
% weight_modi=weight;
% weight_modi(1)=0;
% % weight_modi(2)=0;
% % weight_modi(3)=0;
% weight_modi(N)=0;
% % weight_modi(N-1)=0;
% % weight_modi(N-2)=0;
% % weight_modi(N+1)=0;
% % weight_modi(N+2)=0;
% % weight_modi(end)=0;
% % weight_modi(end-1)=0;
% weight = weight_modi;



Lmax = abs(L*cot(s(1)));

f=f(1:N);
mass_ini = sum(f.*weight')*ds;
counter = 0;
xs=xs(1:N);
f_vec = f;
mass_vec = [0];
t_vec=[0];
RE =[0];
conver=1;
while t<T && conver>1e-5
    mass_error_evo=real(sum(reshape(f,N,1).*weight')*ds-mass_ini)/2;   
	f_tail = f(1:N/2);
	xst=xs(1:N/2);
    options = odeset('RelTol',5e-6,'AbsTol',5e-6);
    fold=f;
    [tt,f_next] = ode15s(@(t,f) odefcn(t,f,P), [0 dt], f,options);
    [f_temp,det,idx] =get_conver_f(f_next,1e-8);
    if det==1
        f=f_temp;
        break
    end
    if norm(f-f_temp)/norm(f)<1e-8
        f=f_temp;
        break
    end
    f=f_temp;
    %mass_error=real(sum(f.*weight)*ds-mass_ini)
    AP_L1 = sum(abs(f-pi^(-0.5)*(1+xs.^2).^(-1)).*weight)*ds;
    f(f<0) = 0;
	t=t+dt
    counter=counter+1;
    f_vec(:,counter)=f;
    % modify f such that mass conserve
    mass_er =real(sum(reshape(f,N,1).*weight')*ds-mass_ini);
    %f(1)=f(1)-mass_er/weight(1)/ds/2; f(end)=f(end)-mass_er/weight(N)/ds/2;
    conver=sum(abs(reshape(fold,N,1)-reshape(f,N,1)).*weight')*ds/2
    mass_vec = [mass_vec mass_error_evo];
    t_vec = [t_vec t];
    %filename = ['LFP_alpha_', num2str_decimal(alpha), '_step_', num2str(counter)];
    %save(filename,'f')
	
end
f_inf=f;
if alpha==1
    correct_tail = pi^(-0.5)*(1+xst.^2).^(-1);
else
    correct_tail =(1+xst).^(-(1+alpha));
end

[L_inf_error,id] = max(abs(f'-pi^(-0.5)*(1+xs(1:N).^2).^(-1)));


figure(1)
plot(xs(1:N),f(1:N),'ro',xs(1:N),pi^(-0.5)*(1+xs(1:N).^2).^(-1),'b','Linewidth',2)
legend('f_{\infty}^{approx}','f_{\infty}^{exact}')
xlim([-5,5])
xlabel('v')
ylabel('f')
title('Equilibrium')
set(gca,'FontSize',35)
set(gcf,'position',[1,1,1440,900])
figurename = ['lf_hom_s_', num2str_decimal(alpha/2), '_final'];
saveas(gcf,figurename,'epsc')
figure(2)
f_tail = f(1:N/2);
xst=xs(1:N/2);
loglog(xst,f_tail,'r-*',xst,correct_tail,'b-^','Linewidth',2)
%loglog(xst,f_tail,'r-*',xst,pi^(-0.5)*(1+xst.^2).^(-1),'b-^')
legend('f_{\infty}','|1+v|^{-(1+2s)}')
%legend('f_{\infty}','f_{\infty}^{exact}')
xlabel('log v')
ylabel('log f')
%title(['Tail shape with T=', num2str(t)])
title('Tail')
set(gca,'FontSize',35)
set(gcf,'position',[1,1,1440,900])
figurename = ['lf_hom_s_', num2str_decimal(alpha/2), '_ts'];
saveas(gcf,figurename,'epsc')
mass_error=real(sum(f.*weight')*ds-mass_ini)/2;

figure(3)
plot(t_vec(2:end),mass_vec(2:end),'Linewidth',2)
%title('Mass error in time')
title(['s=',num2str(alpha/2)])
xlabel('t')
%ylabel('ME')
ylim([-1e-2,1e-2])
set(gca,'FontSize',35)
set(gcf,'position',[1,1,1440,900])
figurename = ['lf_hom_s_', num2str_decimal(alpha/2), '_me'];
saveas(gcf,figurename,'epsc')

% compute relative entropy
RE = zeros(1,length(t_vec)-1);
for i = 1:length(t_vec)-1
    RE(i) = sum(f_vec(:,i).*log(f_vec(:,i)./f_inf').*weight')*ds/2;
end
reft = exp(-2*t_vec(2:end));
figure(4)
%semilogy(t_vec(2:end),RE,'r',t_vec(2:end),reft,'b','Linewidth',2)
semilogy(t_vec(2:end),RE,'r','Linewidth',2)
title('Relative entropy vs time')
%legend('Relative entropy decay rate', 'e^{-2t}')
xlabel('t')
ylabel('log(H^n)')
set(gca,'FontSize',35)
set(gcf,'position',[1,1,1440,900])
figurename = ['lf_hom_s_', num2str_decimal(alpha/2), '_re'];
saveas(gcf,figurename,'epsc')

final_mass_error=mass_vec(end)
filename = ['LF_hom_alpha_', num2str_decimal(alpha),'_N_', num2str(N), '_final'];
save(filename,'f','xs','t','dt','N','weight','mass_error','f_vec','final_mass_error')
end

function dy = odefcn(t,y,P)
dy= P*y;
end



function [P] = LFP_homo_mat(N,alpha,L,l_lim,dt,epsi)
%% Matrix for Levy-fokker-planck operator
s=(pi/N/2:pi/N:2*pi-pi/2/N);
xs=L*cot(s);
ADV = myadv_mat(N);
ADV=ADV(1:N,1:N)+ADV(1:N,N+1:2*N);
FL=fl_mat(N,L,alpha,l_lim);
FL=FL(1:N,1:N)+FL(1:N,N+1:2*N);
A = epsi^(alpha)/dt/2*eye(N);
B=FL-ADV+epsi^(alpha)/dt/2*eye(N);
%P=B^-1*A+eye(N);
%Q=B^(-1)*epsi^(alpha)/dt;
%precond=0;
L=epsi^(-alpha)*(ADV-FL);
%P = eye(N) + dt/2*(L+(eye(N)-dt/2*L)^(-1)*(L+dt/2*L^2));
%P = (eye(N) - dt*epsi^(-alpha)*(ADV-FL))^(-1);
P=epsi^(-alpha)*(ADV-FL);
% M = eye(N) + epsi^(-alpha)*dt*B;
% C=A+B;
% %P=(eye(N)+B^(-1)*A);
% %Q =epsi^(alpha)/dt* B^(-1);
coe1 = ((1+xs(2))/(1+xs(1)))^((1+alpha));
coe2 = ((1+xs(3))/(1+xs(2)))^((1+alpha));
% % M(1,:)=M(2,:)*coe1;
% % M(N,:)=M(N-1,:)*coe1;
% D = eye(N) + epsi^(-alpha/2)*dt^0.5*B;
% %P=B^-1*C;
% %Q = epsi^(alpha)/dt*B^-1;
% precond = D^-1;
% P=M;
% Q=1;
P(2,:)=P(3,:)*coe2;
P(N-1,:)=P(N-2,:)*coe2;
% 
P(1,:)=P(2,:)*coe1;
P(N,:)=P(N-1,:)*coe1;
% % Q(1,:)=Q(2,:)*coe1;
% % Q(N,:)=Q(N-1,:)*coe1;
% % P=precond*M;
% % Q=precond;
end

function [f_temp,det,idx] = get_conver_f(f_next,delta)
[r,c] = size(f_next);
det=0;
for i = 1:r-1
    if norm(f_next(i,:)-f_next(i+1,:))/norm(f_next(i,:))<delta
        f_temp=f_next(i,:)';
        det=1;
        idx=i;
        return
    end
end
f_temp=f_next(end,:);
idx=length(f_next);
end




function FL = fl_mat(N,L,alpha,l_lim)
%filename=['Mk_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*alpha), '_L_', num2str(L)];
%filename=['Mk_pl_alp_1_N_', num2str(N), '_llim_', num2str(l_lim), '_L_', num2str_decimal(L)];
%filename=['Mk_pl_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*alpha), '_L_', num2str_decimal(L)];
%filename=['Mk_plnew_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*alpha), '_L_', num2str_decimal(L)];

if abs(alpha-1)<1e-2
    %filename=['Mk_pl_alp_1_N_', num2str(N), '_llim_', num2str(l_lim), '_L_', num2str_decimal(L)];
    filename=['Mk_plnew_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*alpha), '_L_', num2str_decimal(L)];
else
    filename=['Mk_plnew_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*alpha), '_L_', num2str_decimal(L)];
end



load(filename)
%M=MM(1:N,1:N)+MM(1:N,N+1:2*N);
M=MM;
A=myfft_mat(N);
MA=M*A;
%R=MA(1:N,1:N)+MA(1:N,N+1:2*N);
%FL=(R);
FL=MA;
end


function M=myfft_mat(N)
M=zeros(2*N,2*N);
for j=0:2*N-1
    for k = -N:N-1
        M(j+1,k+1+N) = exp(pi*1i/2/N*k*(2*j+1));
    end
end
M=M^(-1);
end



function ADV=myadv_mat(N)
s=(pi/N/2:pi/N:2*pi-pi/2/N);
%k = 2*[0:N/2-1, 0, -N/2+1:-1];
A=myfft_mat(N);
k = 1*(-N:1:N-1);
K=diag(1i*k);

B=A^(-1);

C=diag(-cos(s).*sin(s));

ADV=(C*B*K*A+eye(2*N));
end


function name=num2str_decimal(a)
s=num2str(a);
c='';
for i = 1:length(s)
    if s(i)=='0'
        c(i)='z';
    elseif s(i)=='.'
        c(i)='p';
    elseif s(i)=='-'
        c(i)='n';
    else
        c(i)=s(i);
    end
end
name=c;
end



