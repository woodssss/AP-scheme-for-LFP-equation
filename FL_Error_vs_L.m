function FL_Error_vs_L(s,N,l_lim,IC)
alpha=2*s;
NN=3;
E = zeros(1,NN); E2=E;
L_vec = [1 3 5 8 10];
%L_vec = [10 12 15];
for ii=1:NN
    L = L_vec(ii);
    compute_fl(N,L,alpha,l_lim,IC);
    %filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    load(filename)
    E(ii) = L_inf_error;
end
alpha=2*alpha;
for ii=1:NN
    L = L_vec(ii);
    compute_fl(N,L,alpha,l_lim,IC);
    %filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    load(filename)
    E2(ii) = L_inf_error;
end
semilogy(L_vec,E,'r',L_vec,E2,'b','Linewidth',2)
title('Error vs Lv')
legend('s=0.2','s=0.4')
xlabel('Lv')
ylabel('E_{\infty}')
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