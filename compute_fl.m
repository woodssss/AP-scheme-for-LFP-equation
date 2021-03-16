function compute_fl(N,L,alpha,l_lim,IC)
close all
ds=pi/N;
%s=(pi/N/2:pi/N:pi*(4*N-1)/2/N);
s=(pi/N/2:pi/N:2*pi-pi/2/N);
xs=L*cot(s);

switch IC
    case 1
        u_0 = exp(-(L*cot(s)).^2);
        fl_exact = 2^alpha*gamma(0.5+alpha/2)/pi^0.5*hypergeom(0.5+alpha/2,0.5,-xs.^2);
        disp('exponential decay IC')
        titlename = ['Exponential decay with s=', num2str(alpha/2)];
    case 2
        u_0=1./(1+(L*cot(s)).^2).^((1-alpha)/2);
        fl_exact = 2^alpha*gamma(0.5+alpha/2)*gamma(0.5-alpha/2)^(-1)*(1+xs.^2).^(-(1+alpha)/2);
        disp('Power law IC with -(1-alpha)')
        titlename = ['Power law decay with s=', num2str(alpha/2)];
    case 3
        u_0=1./(1+(L*cot(s)).^2).^((1+alpha)/2);
        disp('Power law IC with -(1+alpha)')
    otherwise
        error('Case not available, choose from {1,2,3}')
end
u=u_0;
weight=L./((sin(s)).^2);


FTM=myfft_mat(N);
filename=['Mk_plnew_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*alpha), '_L_', num2str_decimal(L)];
load(filename)
M1=MM;
x=(-30:0.2:30);


FL =M1*FTM;
fl_u = FL*u';

figure(1)
plot(xs(1:N),fl_u(1:N),'r-*',xs(1:N),fl_exact(1:N),'b--','Linewidth',2)
legend('u_{approx}', 'u_{exact}')
%xlim([-20,20])
set(gca,'FontSize',35)
xlabel('v')
ylabel('u(v)')
title(titlename)
set(gcf,'Position',[10 10 1500 1000])
L_inf_error = max(abs(fl_u - fl_exact'))
filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
save(filename)


%[L1_error,idx] = max(abs(fl_u - fl_exact'))
%L1_error = sum(abs(fl_u-fl_exact').*weight')*ds
%L2_error = (sum((fl_u - fl_exact').^2)*ds)^0.5
%mass_conserve = sum(fl_u.*weight')*ds
%mass_real = sum(fl_exact.*weight)*ds
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