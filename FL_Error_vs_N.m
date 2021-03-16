function FL_Error_vs_N(s,L,l_lim,IC)
alpha=2*s;
NN=5;
E = zeros(1,NN); E2=E; E3=E;
NL = zeros(1,NN);
for ii=1:NN
    N = 32*2^(ii-1);
    compute_fl(N,L,alpha,l_lim,IC);
    %filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    load(filename)
    E(ii) = L_inf_error;
    NL(ii) = N;
end
alpha=2*alpha;
for ii=1:NN
    N = 32*2^(ii-1);
    compute_fl(N,L,alpha,l_lim,IC);
    %filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    load(filename)
    E2(ii) = L_inf_error;
    NL(ii) = N;
end
alpha=1.6;
for ii=1:NN
    N = 32*2^(ii-1);
    compute_fl(N,L,alpha,l_lim,IC);
    %filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    filename=['compute_fl_alpha_',num2str_decimal(alpha), '_N_', num2str(N), '_L_', num2str_decimal(L), '_IC_', num2str(IC)];
    load(filename)
    E3(ii) = L_inf_error;
    NL(ii) = N;
end
semilogy(NL,E,'r',NL,E2,'b',NL,E3,'c','Linewidth',2)
%semilogy(NL,E,'r',NL,E2,'b','Linewidth',2)
title('Error vs number of mode')
legend('s=0.4','s=0.6','s=0.8')
%legend('s=0.2','s=0.4')
xlabel('Nv')
ylabel('e_{\infty}')
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