function LFP_Energy(Nx,Nv,s,dt,T,epsi,IC)
filename=['LFP_NHE_ssn_alpha_', num2str_decimal(2*s),'_epsi_',num2str_decimal(epsi), '_Nv_',num2str(Nv), '_Nx_', num2str(Nx), '_dt_', num2str_decimal(dt),'_T_',num2str_decimal(T), '_IC_', num2str(IC)];
load(filename)
ll=length(FOM);
xk=(1:1:ll);

plot(xk,GOM.^(0.5),'r-o',xk,ETTM.^(0.5),'b-*',xk,FOM.^(0.5),'c--p','Linewidth',2)
if epsi==1
    temp_s=num2str(1);
else
    temp_s = regexprep(cellstr(num2str(epsi.', '%.0e')), '(?<=e[-+])0*', '');
    temp_s=temp_s{1};
end
title(['Energy stability with s=', num2str(s), ' \epsilon=', temp_s])
legend('E_g', 'E_{\eta}', 'E_f')
xlabel('t')
xk_step = 2*round(ll/10);
label_seq = round((0:xk_step:ll-1)*dt,3);
tick_seq = (1:xk_step:ll);
xticks(tick_seq);
set(gca, 'xticklabel', label_seq);
set(gca,'FontSize',38)
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