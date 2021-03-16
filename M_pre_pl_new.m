function M_pre_pl_new(N,L,alpha,l_lim)
tic

dirname = 'FL_matrix';
if ~exist(dirname, 'dir')
    mkdir(dirname)
end
MM = mk_mat(N,l_lim,alpha,L);
%cd(dirname)
filename=['Mk_plnew_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*alpha), '_L_', num2str_decimal(L)];
save(filename)
%cd ..
time =toc
end




function c=get_c_alpha(a)
c=a*(2^(a-1)*gamma(0.5+a/2))/pi^0.5/gamma(1-a/2);
end


function M = mk_mat(N,l_lim,alpha,big_l)
M = zeros(2*N,2*N);
one_vec = ones(1,2*l_lim+1);
C_alpha=get_c_alpha(alpha);


% compute prod
l_max = 2*N+N*l_lim+N;
prod_seq_common = get_prod_seq_common(alpha, l_max);
prod_seq_even = get_prod_seq_even(alpha, l_max);
prod_seq_odd = get_prod_seq_odd(alpha, l_max);

gn1 = gamma((-1+alpha)/2);
g3 = gamma((3-alpha)/2);

gn = gamma((-alpha)/2);
g2 = gamma((4+alpha)/2);

gn1n = gamma((-1-alpha)/2);
g3p = gamma((3+alpha)/2);



parfor k = 1:2*N %-N:N-1
    temp_v=zeros(2*N,1);
    for j = 1:2*N %0:2*N-1
        sj=pi*(2*(j-1)+1)/2/N;
        if rem(k-N-1,2)==0
            Ak = Ak_even(k-N-1,N,l_lim,alpha,gn1,g3,gn1n,g3p,prod_seq_common,prod_seq_even);
            C = C_alpha*abs(sin(sj))^(alpha-1)/8/big_l^alpha/tan(pi*alpha/2);
        else
            Ak = Ak_odd(k-N-1,N,l_lim,alpha,gn1,g3,gn,g2,prod_seq_common,prod_seq_odd);
            C = 1i*C_alpha*abs(sin(sj))^(alpha-1)/8/big_l^alpha;
        end
        e_vec = get_e_vec(N,j-1);
        temp_v(j) = C*one_vec*Ak*e_vec;
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


function A = Ak_odd(k,N,l_lim,alpha,gn1,g3,gn,g2,prod_seq_common,prod_seq_odd)
A =zeros(2*l_lim+1,N);
for l1 = -l_lim:l_lim
    for l2 = -N/2:N/2-1
        
        
        real_l = abs(l1*N+l2);
        real_lk = abs(k/2-l1*N-l2)-1/2;
        q1 = gamma_quotient_common(gn1,g3,real_l,prod_seq_common);
        q2 = gamma_quotient_odd(gn,g2,real_lk,prod_seq_odd);
        G=q1*q2;
        if isnan(G)
            keyboard;
        end
        A(l1+l_lim+1,l2+N/2+1) = (-1)^l1*((1-alpha)*k^2-4*k*(l1*N+l2))*G*sign(k/2-l1*N-l2);
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