function M_pre_pl_alp1(N,L,l_lim)
tic

dirname = 'FL_matrix';
if ~exist(dirname, 'dir')
    mkdir(dirname)
end
MM = mk_mat(N,l_lim,L);
%cd(dirname)
filename=['Mk_plnew_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*1), '_L_', num2str_decimal(L)];
%filename=['Mk_pl_alp_1_N_', num2str(N), '_llim_', num2str(l_lim), '_L_', num2str_decimal(L)];
save(filename)
%cd ..
time =toc
end




function c=get_c_alpha(a)
c=a*(2^(a-1)*gamma(0.5+a/2))/pi^0.5/gamma(1-a/2);
end


function M = mk_mat(N,l_lim,big_l)
M = zeros(2*N,2*N);
one_vec = ones(1,2*l_lim+1);
%par
for k = 1:2*N %-N:N-1
    temp_v=zeros(2*N,1);
    for j = 1:2*N %0:2*N-1
        sj=pi*(2*(j-1)+1)/2/N;
        if rem(k-N-1,2)==0 % even case
            temp_v(j) = abs(k-N-1)*(sin(sj))^2/big_l*exp(1i*(k-N-1)*sj);
        else
            Bk = Bk_odd(k-N-1,N,l_lim);
            e_vec = get_e_vec(N,j-1);
            temp_v(j) = 1i*(k-N-1)/big_l/pi*(-2/((k-N-1)^2-4)-one_vec*Bk*e_vec);
        end
        %temp_v(j) = C*one_vec*Ak*e_vec;
        %temp_v(j)
    end
    v=reshape(temp_v,2*N,1);
    M(:,k)=v;
end
end


function e_vec = get_e_vec(N,j)
e_vec = zeros(N,1);
sj=pi*(2*j+1)/2/N;
for l = -N/2:N/2-1
    e_vec(l+N/2+1) = exp(1i*2*l*sj);
end

end

function A = Ak_even(k,N,l_lim,alpha,gn1,g3,gn1n,g3p,prod_seq_common,prod_seq_even)
A =zeros(2*l_lim+1,N);
for l1 = -l_lim:l_lim
    for l2 = -N/2:N/2-1
        
        
        real_l = abs(l1*N+l2);
        real_lk = abs(k/2-l1*N-l2);
        q1 = gamma_quotient_common(gn1,g3,real_l,prod_seq_common);
        q2 = gamma_quotient_even(gn1n,g3p,real_lk,prod_seq_even);
        
        G = q1*q2;    
        if isnan(G)
            keyboard;
        end
        A(l1+l_lim+1,l2+N/2+1) = (-1)^l1*((1-alpha)*k^2-4*k*(l1*N+l2))*G;
    end
end

end


function A = Bk_odd(k,N,l_lim)
A =zeros(2*l_lim+1,N);
for l1 = -l_lim:l_lim
    for l2 = -N/2:N/2-1
        
        
        top = 4*(-1)^l1*sign(l1*N+l2);
        bot = (k-2*(l1*N+l2))*((k-2*(l1*N+l2))^2-4);
        G=top/bot;
        if isnan(G)
            keyboard;
        end
        A(l1+l_lim+1,l2+N/2+1) = G;
    end
end

end

% functions for common terms

function g = gamma_quotient_common(gn1,g3,real_l,prod_seq_com)
% gn1 = gamma((-1+alpha)/2)
% g3 = gamma((3-alpha)/2)
% prod_seq is prod_{m=0}^{|l|-1}

if real_l ==0
    g= gn1/g3;
else
    g= gn1/g3*prod_seq_com(real_l);
end
end

function prod_seq = get_prod_seq_common(alpha, l_max)
seq = zeros(1,l_max);
for m = 0:l_max-1
    seq(m+1) = ((-1+alpha)/2 + m)/((3-alpha)/2+m);
end
prod_seq = seq;
for i = 2:l_max
    prod_seq(i) = prod_seq(i)*prod_seq(i-1);
end
end



% functions for even terms

function g = gamma_quotient_even(gn1n,g3p,real_l,prod_seq_even)
% gn1n = gamma((-1-alpha)/2)
% g3p = gamma((3+alpha)/2)
% prod_seq is prod_{m=0}^{|l|-1}

if real_l == 0
    g=gn1n/g3p;
else
    g= gn1n/g3p*prod_seq_even(real_l);
end
end

function prod_seq = get_prod_seq_even(alpha, l_max)
seq = zeros(1,l_max);
for m = 0:l_max-1
    seq(m+1) = ((-1-alpha)/2 + m)/((3+alpha)/2+m);
end
prod_seq = seq;
for i = 2:l_max
    prod_seq(i) = prod_seq(i)*prod_seq(i-1);
end
end


% functions for odd terms

function g = gamma_quotient_odd(gn,g2,real_l,prod_seq_odd)
% gn = gamma((-alpha)/2)
% g2 = gamma((4+alpha)/2)
% prod_seq is prod_{m=0}^{|l|-1}
if real_l ==0
    g= gn/g2;
else
    g= gn/g2*prod_seq_odd(real_l);
end
end

function prod_seq = get_prod_seq_odd(alpha, l_max)
seq = zeros(1,l_max);
for m = 0:l_max-1
    seq(m+1) = ((-alpha)/2 + m)/((4+alpha)/2+m);
end
prod_seq = seq;
for i = 2:l_max
    prod_seq(i) = prod_seq(i)*prod_seq(i-1);
end
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