function LFP_FOT_check(Nx,Nv,Lx,Lv,s,l_lim,gamma,epsi,T,IC)
%%%%%%%%%%%%% Check AP for 1D case %%%%%%
% This solver aims to solve
% \epsi^2s \patial_t f + \epsi v \nabla_x f = LFP(f).
% tau is the preset outer time step
% Nx Nv are the number of uniform grid point
% T is final time
% Author: Wuzhe Xu
% Date: 10/15/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NNN=5;
legend_name = cell(1,NNN);
if s<0.5
    st=2;
else
    st=2;
end
rho_t_vec = zeros(NNN,Nx);
for ii = 1:NNN
    dt = T*2^(1-ii);
    %LFP_AP_HE(Nx,Nv,Lx,Lv,s,epsi,l_lim,dt,T,2);
    LFP_AP_NHE_ssn(Nx,Nv,Lx,Lv,s,epsi,l_lim,gamma,dt,T,IC,0);
    %PN_search_VPFP_1D(temp_epsi,Nx,Nv,Lx,Lv,tau,T);
    %filename=['LFP_alpha_', num2str_decimal(2*s),'_epsi_',num2str_decimal(epsi),'_dt_', num2str_decimal(dt),'_T_',num2str_decimal(T)];
    %filename=['PN_vpfp_epsi', num2str_decimal(temp_epsi),'_t_'   ,  num2str_decimal(T)];
    filename=['LFP_NHE_ssn_alpha_', num2str_decimal(2*s),'_epsi_',num2str_decimal(epsi), '_Nv_',num2str(Nv), '_Nx_', num2str(Nx), '_dt_', num2str_decimal(dt),'_T_',num2str_decimal(T), '_IC_', num2str(IC)];
    load(filename)
    rho_t_vec(ii,:) = rho_approx;
%     E2M_vec(ii,:) = E2M;
%     temp_s = regexprep(cellstr(num2str(epsi.', '%.0e')), '(?<=e[-+])0*', '');
%     legend_name{ii}=['\epsilon =', temp_s{1}];
%     ii
end
ET=zeros(1,NNN-1);
for jj = 1:NNN-1
    %ET(jj) = (sum((abs(real(rho_t_vec(jj,:))-real(rho_t_vec(jj+1,:)))).^2))^0.5*dx;
    ET(jj) = (sum((abs(real(rho_t_vec(jj,:))-real(rho_t_vec(jj+1,:))))))*dx;
    %ET(jj) = norm(rho_t_vec(jj,:)-rho_t_vec(jj+1,:));
    %(sum((abs(real(rho_t_vec(jj,:))-real(rho_t_vec(jj+1,:)))).^2))*dx;
end
xk=flip([T/2 T/4 T/8 T/16]);
%xk=flip([0.1 0.05 ]);
loglog(xk,flip(ET),'r-*',xk,xk,'b--')
% semilogy(xk,E2M_vec(1,1:end-1),'r-o',xk,E2M_vec(2,1:end-1),'b-^',xk,E2M_vec(3,1:end-1),'g-x',xk,E2M_vec(4,1:end-1),'m-p','Linewidth',2)
title('1D LFP: 1st order accuracy in time')
%legend(legend_name)
xlabel('\Delta t','Fontsize',25)
ylabel('e_{\Delta t}','Fontsize',20)
set(gca,'FontSize',35)
set(gcf,'position',[1,1,1440,900])
filename = ['FOT_alpha_',num2str_decimal(2*s),'_epsi_',num2str_decimal(epsi),'_dt_', num2str_decimal(dt),'_Nx_',num2str(Nx),'_Nv_',num2str(Nv),'_T_',num2str_decimal(T),'_IC_',num2str(IC)];
save(filename,'xk','T','ET')
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