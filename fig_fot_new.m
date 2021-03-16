function fig_fot_new(Nx,Nv,s,T,IC)
dt=T/2^4/4;
epsi1=1;
filename = ['FOT_f_alpha_sdt_',num2str_decimal(2*s),'_epsi_',num2str_decimal(epsi1),'_dt_', num2str_decimal(dt),'_Nx_',num2str(Nx),'_Nv_',num2str(Nv),'_T_',num2str_decimal(T),'_IC_',num2str(IC)];
load(filename)
ET1 = ET;

epsi3=1e-3;
filename = ['FOT_f_alpha_sdt_',num2str_decimal(2*s),'_epsi_',num2str_decimal(epsi3),'_dt_', num2str_decimal(dt),'_Nx_',num2str(Nx),'_Nv_',num2str(Nv),'_T_',num2str_decimal(T),'_IC_',num2str(IC)];
load(filename)
ET3 = ET;

epsi5=1e-5;
filename = ['FOT_f_alpha_sdt_',num2str_decimal(2*s),'_epsi_',num2str_decimal(epsi5),'_dt_', num2str_decimal(dt),'_Nx_',num2str(Nx),'_Nv_',num2str(Nv),'_T_',num2str_decimal(T),'_IC_',num2str(IC)];
load(filename)
ET5 = ET;
xk = flip([0.025 0.0125 0.00625 0.003125]);
loglog(xk,flip(ET1),'b-*',xk,flip(ET3),'k-p',xk,flip(ET5),'m-o',xk,xk,'r--','Linewidth',2)
% semilogy(xk,E2M_vec(1,1:end-1),'r-o',xk,E2M_vec(2,1:end-1),'b-^',xk,E2M_vec(3,1:end-1),'g-x',xk,E2M_vec(4,1:end-1),'m-p','Linewidth',2)
%title('1D VLFP: 1st order accuracy in time')
legend('\epsilon = 1','\epsilon = 1e-3','\epsilon = 1e-5', 'slope = 1','Location','northwest')
xlabel('\Delta t','Fontsize',25)
ylabel('e_{\Delta t}','Fontsize',20)
xlim([0.0025, 0.03])
ylim([1e-4, 1e-1])
set(gca,'FontSize',35)
set(gcf,'position',[1,1,1440,900])
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