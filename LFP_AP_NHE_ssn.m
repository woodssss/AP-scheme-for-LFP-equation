function LFP_AP_NHE_ssn(Nx,Nv,Lx,Lv,s,epsi,l_lim,gamma,dt,T,IC,det1)
%%%%%%%%%%%%%%%%%% Asymptotic preserving scheme for LFP %%%%%%%%%%%%
% Asymptotic preserving scheme with Novel Hilbert expansion
% consider \epsi^(2s) \partial_t f + \epsi v \cdot \nabla_x f
%                                   = \nabla \cdot (vf)- (-\Delta^s_v) f.
%
% We use change of variable on v direction, because of fat tail
% equilibrium.
%
% Use [-Lx,Lx] as domain and partition into Nx pieces uniformly.
%
% Nv Lv l_lim are coefficients in change of variable and computation of FL.
%
% alpha = 2s.
% 
% dt is the time step and T is the final time.
%
% This scheme preserves mass and gives correct limit system.
%
% On x direction, we use periodic BC, thus the FL on x direction can be
% computed by FFT.
%
% Consider f0 = pi^(-0.5)*(1+sin(2*pi*(x/Lx)))*exp(-v^2).
%
% Wuzhe Xu, 10/14/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%% define coefficient
% Set IC
if nargin<11
    det1 = 0;
end
alpha=2*s;
dx = 2*Lx/(Nx);
x = (-Lx+dx/2:dx:Lx-dx/2); xx=(-Lx+dx/4:dx/2:Lx-dx/4);

ds = pi/Nv;
ss=(pi/Nv/2:pi/Nv:2*pi-pi/2/Nv)';
dds=pi/2/Nv; sss=(pi/Nv/4:pi/Nv/2:2*pi-pi/2/Nv/2)'; vvs=Lv*cot(ss(1:2*Nv));
vs=Lv*cot(ss(1:Nv));
Ns=Nv;
weight=Lv./((sin(ss)).^2);
weight=weight(1:Nv);
weight(1)=0;weight(end)=0;

%Define equilibrium 

% get M_real
filename1=['LF_hom_alpha_',num2str_decimal(2*s),'_N_', num2str(Nv),'_final'];
load(filename1,'f');
M_real=f(1:Nv)/(sum(f.*weight')*ds);

M = M_real;

% M_temp = 1./(1+abs(vs).^(1+2*s));
% M = M_temp'/(sum(M_temp.*weight)*ds);
% 
% M_temp1 = 1./(1+abs(vs).^(1+2*s));
% M1 = M_temp1'/(sum(M_temp1.*weight)*ds);

%% Initial condition
%rho0 = (1+sin(2*pi*(x/Lx)));
rt0=zeros(Nx,Nv);
f0=rt0;
N=Nx*Nv;
switch IC
    case 1
        a=15;
        rho0 = exp(-a*x.^2);
        for i=1:Nx
            for j=1:Nv
                rt0(i,j) = exp(-a*(x(i)+epsi*vs(j))^2);
                f0(i,j) = (pi)^(-0.5)*exp(-a*x(i)^2)*exp(-vs(j)^2);
            end
        end
        if det1==1
            f0_imx=zeros(2*Nx,Nv);
            for i=1:2*Nx
                for j=1:Nv
                    f0_imx(i,j) = (pi)^(-0.5)*exp(-a*xx(i)^2)*exp(-vs(j)^2);
                end
            end
        end
        if det1==3
            f0_imx=zeros(Nx,Nv);
            for i=1:Nx
                for j=1:Nv
                    f0_imx(i,j) = (pi)^(-0.5)*exp(-a*x(i)^2)*exp(-vs(j)^2);
                end
            end
        end
        rho_real = exp(-a*x.^2);
        disp('exponential decay IC')
    case 2
        rho0 = (1+sin(pi*(x/Lx)));
        for i=1:Nx
            for j=1:Nv
                rt0(i,j) = (1+sin(1*pi*((x(i)+epsi*vs(j))/Lx)));
                f0(i,j) =pi^(-0.5)*(1+sin(1*pi*(x(i)/Lx)))*exp(-vs(j)^2);
            end
        end
        if det1==1
            f0_imx=zeros(2*Nx,Nv);
            for i=1:2*Nx
                for j=1:Nv
                    f0_imx(i,j) =pi^(-0.5)*(1+sin(1*pi*(xx(i)/Lx)))*exp(-vs(j)^2);
                end
            end
        end
        if det1==3
            f0_imx=zeros(Nx,Nv);
            for i=1:Nx
                for j=1:Nv
                    f0_imx(i,j) =pi^(-0.5)*(1+sin(1*pi*(x(i)/Lx)))*exp(-vs(j)^2);
                end
            end
        end
        rho_real = (1+sin(1*pi*(x/Lx)));
        disp('sin IC')
end
g0 = f0-rt0.*(repmat(M,Nx,1));

%% determine dt
%% Exact rho evo
%% Get Evolution matrix
[FL,ADV] = get_M(Nv,alpha,Lv,l_lim);
DFT = dftmtx(Nx);
k = pi/Lx*[0:Nx/2-1, 0, -Nx/2+1:-1];
K=diag(1i*k);
Ks=diag(abs(k).^(2*s));
IFT = conj(dftmtx(Nx))/Nx;
TEV1 = ((epsi^(2*s)/dt+gamma)*eye(Nv)-ADV+FL)^(-1);
FLEVx = (eye(Nx)+dt*IFT*Ks*DFT)^(-1);
TEV2_cell=cell(1,Nv);
for j=1:Nv
    TEV2_cell{1,j}=((epsi^(2*s)/dt-gamma)*eye(Nx)+epsi*vs(j)*IFT*K*DFT)^(-1);
end
if det1 == 1
%% determine if we want to run IMEX scheme
    kk = pi/Lx*[0:Nx-1, 0, -Nx+1:-1];
    K2=diag(1i*kk);
    DFT2 = dftmtx(2*Nx);
    IFT2=conj(dftmtx(2*Nx))/Nx/2;
    ADXimx=IFT2*K2*DFT2;
    sub_dt=dt/50;
    EVimx = (eye(Nv)-epsi^(-2*s)*sub_dt*(ADV-FL))^(-1);
    %EVimx = (eye(Nv)+epsi^(-2*s)*sub_dt*(ADV-FL));
    ts_cell=cell(1,Nv);
    for j=1:Nv
        ts_cell{1,j}=(eye(2*Nx)+epsi^(1-2*s)*vs(j)*sub_dt*ADXimx)^(-1);
        %ts_cell{1,j}=(eye(2*Nx)-epsi^(1-2*s)*vs(j)*sub_dt*ADXimx);
    end
end

if det1 == 3
%% determine if we want to run IMEX scheme
    kk = pi/Lx*[0:Nx/2-1, 0, -Nx/2+1:-1];
    K2=diag(1i*kk);
    DFT2 = dftmtx(Nx);
    IFT2=conj(dftmtx(Nx))/Nx;
    ADXimx=IFT2*K2*DFT2;
    sub_dt=dt/50;
    EVimx = (eye(Nv)-epsi^(-2*s)*sub_dt*(ADV-FL))^(-1);
    %EVimx = (eye(Nv)+epsi^(-2*s)*sub_dt*(ADV-FL));
    ts_cell=cell(1,Nv);
    for j=1:Nv
        ts_cell{1,j}=(eye(Nx)+epsi^(1-2*s)*vs(j)*sub_dt*ADXimx)^(-1);
        %ts_cell{1,j}=(eye(Nx)-epsi^(1-2*s)*vs(j)*sub_dt*ADXimx);
    end
end

t=0;
g=g0;
trho=rt0;
f_vec=reshape(f0',N,1);
ct=2;
E2M = zeros(1,T/dt+1); E2imx=E2M;
f_limit = repmat(rho_real',1,Nv).*repmat(M,Nx,1);
E2M(1,1)=sum(sum(abs(real(   (f0-f_limit).*repmat(weight',Nx,1) ))))*ds*dx;
E2imx(1,1)=0;
if det1==1 || det1==3
    f_imx=f0_imx;
end
GOM=zeros(1,T/dt+1); ETTM=zeros(1,T/dt+1); FOM=zeros(1,T/dt+1);;
GOM(1,1)=sum(sum(abs(real((g0.^2./(repmat(M,Nx,1))).*repmat(weight',Nx,1) ))))*ds*dx;
ETTM(1,1) = sum(sum(abs(real((trho.^2.*(repmat(M,Nx,1))).*repmat(weight',Nx,1) ))))*ds*dx;
FOM(1,1) = sum(sum(abs(real((f0.^2./(repmat(M,Nx,1))).*repmat(weight',Nx,1) ))))*ds*dx;
while abs(t-T)>1e-3
    %% macro-micro decomposition
    % evo of g
    if det1==3
        g_star = evo_g_1ss(trho,g,M,ADV,FL,TEV1,Lx,Nx,Nv,s,epsi,dt);
    else
        g_star = evo_g_1ss(trho,g,M,ADV,FL,TEV1,Lx,Nx,Nv,s,epsi,dt);
    end
    gn=evo_g_2ss(g_star,Nv,s,epsi,dt,TEV2_cell);
    %norm(gn)
    trhon = evo_trho_im_ss(trho,FLEVx,Nv);
    %trhon = evo_trho_ex(trho,Lx,Nx,Nv,s,dt);
    trho=trhon;
    g=gn;
    
    % recover f
    fmm=trhon.*(repmat(M,Nx,1))+gn;
    %fmm=trhon.*(repmat(M,Nx,1));
    fmm(fmm<1e-10)=1e-10;
    
    % calculate rho_approx
    rho_approx = average(fmm,Nx,weight,ds);
    
    %% IMEX
    if det1 ==1
        t_temp=0; 
        while t_temp<dt
            % transport
            for j=1:Nv
                f_imx(:,j)=ts_cell{1,j}*f_imx(:,j);
                %f_imx(:,j)=(eye(Nx)-epsi^(1-2*s)*vs(j)*sub_dt*ADXimx)*f_imx(:,j);
            end
            % collision
            for i=1:2*Nx
                f_imx(i,:)=(EVimx*f_imx(i,:)')';
            end
            t_temp=t_temp+sub_dt;
        end
        rho_imx = average(f_imx,2*Nx,weight,ds);
    end
    
    if det1 ==3
        t_temp=0; 
        while t_temp<dt
            % transport
            for j=1:Nv
                f_imx(:,j)=ts_cell{1,j}*f_imx(:,j);
                %f_imx(:,j)=(eye(Nx)-epsi^(1-2*s)*vs(j)*sub_dt*ADXimx)*f_imx(:,j);
            end
            % collision
            for i=1:Nx
                f_imx(i,:)=(EVimx*f_imx(i,:)')';
            end
            t_temp=t_temp+sub_dt;
        end
        rho_imx = average(f_imx,Nx,weight,ds);
    end
    
    %% Calculate limit rho
    rho_real = rho_evo_new_im(rho_real,dt,Nx,Lx,alpha);
    %rho_real = rho_evo_new_ex(rho_real,dt,Nx,Lx,s);
    f_limit = repmat(rho_approx',1,Nv).*repmat(M_real,Nx,1);
    mass_rho = sum(rho_approx)*dx;
    mass_rho_real = sum(rho_real)*dx;
    %error2real = sum(abs(real(rho_approx)-rho_real))*dx;
    %error2real=sum(sum(abs(real(   (fmm-f_limit).*repmat(weight',Nx,1) ))))*ds*dx;
    error2real=sum(sum(abs(real((fmm-f_limit).*repmat(weight',Nx,1) ))))*ds*dx;
    %e1= sum(abs(real(rho_approx)-rho_real))*dx
    %e2=  max(max(abs(real(  repmat(real(rho_approx')-rho_real',1,Nv).*repmat(M,Nx,1).*repmat(weight',Nx,1)))))*ds*dx
    %norm(fmm-f_limit)
    E2M(1,ct) = error2real;
    GOM(1,ct) = sum(sum(abs(real((gn.^2./(repmat(M,Nx,1))).*repmat(weight',Nx,1) ))))*ds*dx;
    ETTM(1,ct) = sum(sum(abs(real((trho.^2.*(repmat(M,Nx,1))).*repmat(weight',Nx,1) ))))*ds*dx;
    FOM(1,ct) = sum(sum(abs(real((fmm.^2./(repmat(M,Nx,1))).*repmat(weight',Nx,1) ))))*ds*dx;
    disp([' Mass of rho exact at current step is ', num2str(mass_rho)])
    disp([' Mass of rho approx at current step is ', num2str(mass_rho_real)])  
    if det1==3
        tmp_error=sum(sum(abs(real((fmm-f_imx).*repmat(weight',Nx,1) ))))*ds*dx;
        disp([' L1 error between f and f_imx ', num2str(tmp_error)]) 
        E2imx(1,ct)=tmp_error;
    end
    %tmp_error=sum(sum(abs(real((fmm-f_imx).*repmat(weight',Nx,1) ))))*ds*dx;
     
    t=t+dt;
    ct=ct+1;
end
if det1 ==1
    disp('running two schemes')
    figure(1)
    %plot(x,rho_real,'r--',x,real(rho_approx),'b*',x,rho_ip,'c^')
    %plot(x,real(rho_approx),'b--',x,rho_ip,'r*',x,rho0,'c','Linewidth',2)
    plot(x,real(rho_approx),'b-*',xx,rho_imx,'r--','Linewidth',2)
    drawnow
    %legend('\rho approx by macro-micro scheme', '\rho approx by fully implicit scheme','initial condition','Location','northwest')
    legend('AP scheme', 'reference \rho by IMEX','Location','northwest')
    title(['\epsilon = ', num2str(epsi), ' s = ', num2str(s), ' t=', num2str(t)])
    xlim([-3,3])
    xlabel('x')
    ylabel('\rho')
    set(gca,'FontSize',35)
    set(gcf,'position',[1,1,1440,900])
elseif det1 ==2
    disp('running for AP')
    plot(x,rho_real,'r*',x,real(rho_approx),'b--','Linewidth',2)
    drawnow
    legend('\rho exact', '\rho approx by mirco-macro scheme','Location','northwest')
    title(['\epsilon = ', num2str(epsi), ' t=', num2str(t)])
    xlabel('x')
    ylabel('\rho')
    set(gca,'FontSize',35)
    set(gcf,'position',[1,1,1440,900])
elseif det1==3
    disp('running two schemes')
    figure(1)
    %plot(x,rho_real,'r--',x,real(rho_approx),'b*',x,rho_ip,'c^')
    %plot(x,real(rho_approx),'b--',x,rho_ip,'r*',x,rho0,'c','Linewidth',2)
    plot(x,real(rho_approx),'b-*',x,rho_imx,'r-o','Linewidth',2)
    drawnow
    %legend('\rho approx by macro-micro scheme', '\rho approx by fully implicit scheme','initial condition','Location','northwest')
    legend('\rho by mirco-macro scheme', 'reference \rho by IMEX','Location','northwest')
    title(['\epsilon = ', num2str(epsi), ' t=', num2str(t)])
    %xlim([-3,3])
    xlabel('x')
    ylabel('\rho')
    set(gca,'FontSize',35)
    set(gcf,'position',[1,1,1440,900])
end
filename=['LFP_NHE_ssn_alpha_', num2str_decimal(2*s),'_epsi_',num2str_decimal(epsi), '_Nv_',num2str(Nv), '_Nx_', num2str(Nx), '_dt_', num2str_decimal(dt),'_T_',num2str_decimal(T), '_IC_', num2str(IC)];
save(filename)

end

function I = ITRM(trho,M,P,Nx,Nv,epsi,s,Lx)
% P is LF
I = zeros(Nx,Nv);
k = pi/Lx*[0:Nx/2-1, 0, -Nx/2+1:-1]';
flvtrho = zeros(Nx,Nv);
for j = 1:Nv
    flvtrho(:,j) = ifft(abs(k).^(2*s).*fft(trho(:,j)));
end
flvtrho = epsi^(2*s)*flvtrho;
for i = 1:Nx
   trhom = trho(i,:).*M;
   I(i,:) =  (P*trhom'-(flvtrho(i,:)').*M'-trho(i,:)'.*(P*M'))';
end
end

function gn=evo_g_1ss(trho,g,M,ADV,FL,TEV1,Lx,Nx,Nv,s,epsi,dt)
I = ITRM(trho,M,FL,Nx,Nv,epsi,s,Lx);
%II2 = I - trho.*(repmat(((ADV-FL)*M')',Nx,1));
for i=1:Nx
    %g(i,:) = epsi^(2*s)/dt*(TEV1*g(i,:)')'-(TEV1*II2(i,:)')';
    g(i,:) = epsi^(2*s)/dt*(TEV1*g(i,:)')'-(TEV1*I(i,:)')';
end
gn=g;
end

function gn=evo_g_1ss3(trho,g,M,ADV,FL,TEV1,Lx,Nx,Nv,s,epsi,dt)
I = ITRM(trho,M,FL,Nx,Nv,epsi,s,Lx);
II2 = I - trho.*(repmat(((ADV-FL)*M')',Nx,1));
for i=1:Nx
    g(i,:) = epsi^(2*s)/dt*(TEV1*g(i,:)')'-(TEV1*II2(i,:)')';
    %g(i,:) = epsi^(2*s)/dt*(TEV1*g(i,:)')'-(TEV1*I(i,:)')';
end
gn=g;
end



function gn=evo_g_2ss(g,Nv,s,epsi,dt,TEV2_cell)
for j=1:Nv
    TEV=TEV2_cell{1,j};
    g(:,j)=epsi^(2*s)/dt*TEV*g(:,j);
end
gn=g;
end

function g_vec = m2vecx(g,Nx,Nv)
g_vec = reshape(g,Nx*Nv,1);
end

function g = vecx2m(g_vec,Nx,Nv)
g = reshape(g_vec,Nx,Nv);
end

function trhon = evo_trho_ex(trho,Lx,Nx,Nv,s,dt)
k = pi/Lx*[0:Nx/2-1, 0, -Nx/2+1:-1]';
flxtrho = trho;
for j = 1:Nv
    flxtrho(:,j) = ifft((abs(k)).^(2*s).*fft(trho(:,j)));
end
trhon = trho - dt*flxtrho;
end

function trhon = evo_trho_im_ss(trho,FLEVx,Nv)
for j=1:Nv
    trhon(:,j) = FLEVx*trho(:,j);
end
end

function trhon = evo_trho_im(trho,TREV,Nx,Nv)
trho_vec = m2vecx(trho,Nx,Nv);
trho_vec_n = TREV*trho_vec;
trhon=vecx2m(trho_vec_n,Nx,Nv);
end




function avg = average(f,Nx,weight,ds)
avg = zeros(1,Nx);
for i = 1:Nx
    avg(i) = sum(f(i,:).*weight')*ds;
end
end

function nrho = rho_evo_new_im(rho,dt,Nx,Lx,alpha)
%% Implicit scheme
k = pi/Lx*[0:Nx/2-1, 0, -Nx/2+1:-1];
FT = dftmtx(Nx);
IFT=conj(dftmtx(Nx))/Nx;
Ks = diag((abs(k)).^(alpha));
EV = (eye(Nx)+dt*IFT*Ks*FT)^(-1);
%nrho = ifft((1./(1+dt*abs(k).^alpha)).*fft(rho));
nrho = (EV*rho')';
end

function nrho = rho_evo_new_ex(rho,dt,Nx,Lx,s)
%% Explicit scheme
k = pi/Lx*[0:Nx/2-1, 0, -Nx/2+1:-1];
FT = dftmtx(Nx);
IFT=conj(dftmtx(Nx))/Nx;
Ks = diag((abs(k)).^(2*s));
TEV = (eye(Nx)-dt*IFT*Ks*FT);
nrho=(TEV*rho')';
end

function [FL,ADV] = get_M(N,alpha,L,l_lim)
ADV = new_ad(N);
ADV=ADV(1:N,1:N)+ADV(1:N,N+1:2*N);
FL=fl_mat(N,L,alpha,l_lim);
FL=FL(1:N,1:N)+FL(1:N,N+1:2*N);
end

function AD = new_ad(N)
s=(pi/N/2:pi/N:2*pi-pi/2/N);
%k = 2*[0:N/2-1, 0, -N/2+1:-1];
A=myfft_mat(N);
k = 1*(-N:1:N-1);
%k = [0:N-1, 0, -N+1:-1];
K=diag(1i*k);

B=A^(-1);

C=diag(-cos(s).*sin(s));

AD=(C*B*K*A+eye(2*N));
%AD=(C*B*K*A);

end

function FL = fl_mat(N,L,alpha,l_lim)
%filename=['Mk_N_', num2str(N), '_llim_', num2str(l_lim), '_alpha_', num2str(10*alpha), '_L_', num2str(L)];
if abs(alpha-10)<1e-2
    filename=['Mk_pl_alp_1_N_', num2str(N), '_llim_', num2str(l_lim), '_L_', num2str_decimal(L)];
%elseif abs(alpha-0.4)<1e-2
    %filename=['Mk_N_', num2str(N), '_llim_', num2str(l_lim),'_alpha_4', '_L_', num2str_decimal(L)];
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

function [P,Q] = LFP_homo_mat(N,alpha,L,l_lim,dt,epsi)
%% Matrix for Levy-fokker-planck operator
s=(pi/N/2:pi/N:2*pi-pi/2/N);
xs=L*cot(s);
ADV = myadv_mat(N);
ADV=ADV(1:N,1:N)+ADV(1:N,N+1:2*N);
FL=fl_mat(N,L,alpha,l_lim);
FLO=FL;
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
Q = FLO;
% % Q(1,:)=Q(2,:)*coe1;
% % Q(N,:)=Q(N-1,:)*coe1;
% % P=precond*M;
% % Q=precond;
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