function LFP_AP_check_new(Nx,Nv,Lx,Lv,gamma,s,l_lim,dt,T,IC)
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
    st=1;
else
    st=1;
end
for ii = 1:NNN
    temp_epsi = 10^(-st-ii);
    %LFP_AP_HE(Nx,Nv,Lx,Lv,s,temp_epsi,l_lim,dt,T,1);
    %LFP_AP_NHE(Nx,Nv,Lx,Lv,s,temp_epsi,l_lim,dt,T,IC,2);
    LFP_AP_NHE_ssn(Nx,Nv,Lx,Lv,s,temp_epsi,l_lim,gamma,dt,T,IC,0);
    %PN_search_VPFP_1D(temp_epsi,Nx,Nv,Lx,Lv,tau,T);
    filename=['LFP_NHE_ssn_alpha_', num2str_decimal(2*s),'_epsi_',num2str_decimal(temp_epsi), '_Nv_',num2str(Nv), '_Nx_', num2str(Nx), '_dt_', num2str_decimal(dt),'_T_',num2str_decimal(T), '_IC_', num2str(IC)];
    load(filename)
    E2M_vec(ii,:) = E2M;
    temp_s = regexprep(cellstr(num2str(temp_epsi.', '%.0e')), '(?<=e[-+])0*', '');
    legend_name{ii}=['\epsilon =', temp_s{1}];
    ii
end
ll=length(E2M_vec);
xk=(1:1:ll);
semilogy(xk,E2M_vec(1,1:end),'r-o',xk,E2M_vec(2,1:end),'b-^',xk,E2M_vec(3,1:end),'g-x',xk,E2M_vec(4,1:end),'m-p',xk,E2M_vec(5,1:end),'k-*','Linewidth',2)
title(['1D LFP: asymptotic behavior with s=',num2str(s)])
legend(legend_name)
xlabel('Time','Fontsize',25)
%ylabel('log||f^{\epsilon}-\rho M_{ref}||_1','Fontsize',20)
%ylim([1e-10,1e1])
xk_step = round(ll/10);
label_seq = round((0:xk_step:ll-1)*dt,3);
tick_seq = (1:xk_step:ll);
xticks(tick_seq);
set(gca, 'xticklabel', label_seq);
set(gca,'FontSize',30)
set(gcf,'position',[1,1,1440,900])
filename=['LFP_AP_check_IC', num2str(IC)];
save(filename, 'xk', 'E2M_vec')
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