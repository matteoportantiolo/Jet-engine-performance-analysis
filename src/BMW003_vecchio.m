function [THRUST, TSFC, T0, P0, rho0, T0_t, P0_t, rho0_t, A_in, T1, P1, rho1, a1, v_in, M_volo, M_in, M_e, ...
                v_e, f, m_a, m_f, m_tot, P_p, P_dis, P_av, eta_p, eta_ter, eta_glob, P5_st, m_a2, m_norm,I_sp]...
                = BMW003_vecchio (z, v_volo, pi_d, beta_c, eta_c_tot, eta_c, eta_b, pi_b, deltaH, eta_t_t, eta_t_m, eta_n)

%aria
gamma_a=1.4;
R_a = 287;
cp_a=1004;
%gas combusti
gamma_gc=1.33;
R_gc=286.6;
cp_gc=1155;

%quota
p_amb0 = 101325;
T_amb0 = 288.15;
rho_amb0 = 1.225;
a=0.0065;
g = 9.81;
T0=T_amb0-a*z;
P0=p_amb0*(1-(a*z)/T_amb0).^(g/(R_a*a));
rho0=P0/(R_a*T0);

THRUST = [];
m_a = [];
v_in = [];
v_e = [];
M_e = [];
m_norm = [];
f = [];
eta_p = [];
eta_ter = [];
eta_glob = [];
TSFC = [];
I_sp = [];
M_in = [];

for i = 1:length(v_volo)
v_vol = v_volo(i).*1.609/3.6;
a_volo = sqrt(gamma_a*R_a*T0);
M_volo = v_vol./a_volo;
T0_t = T0.*(1+(gamma_a-1)./2.*M_volo.^2);
P0_t = P0.*(1+(gamma_a-1)./2.*M_volo.^2).^(gamma_a/(gamma_a-1));
rho0_t = rho0.*(1+(gamma_a-1)./2.*M_volo.^2).^(1/(gamma_a-1));
T3_t = 1026.95; %T in camera di combustione

%eta_d = ((1+(gamma_a-1)/2*M_volo^2)*pi_d^((gamma_a-1)/gamma_a)-1)/((gamma_a-1)/2*M_volo^2);    %non l'ho usato; 
% serve se si considera il flusso come non adiabatico (a velocità basse si
% abbassa molto (guardare slide 14 della presa d'aria)
eta_c_m=eta_c_tot/eta_c;

eta_t = eta_t_t/eta_t_m;

%% Area di ingresso
D_in = 20;
D_in = D_in*2.54/100;
A_in = pi*D_in^2/4;

A_e = 0.1020;      %efflusso, la considero costante tanto varia al massimo tra 970 e 1020 cm^2

%% 0--->1: diffusore
P1_t = pi_d.*P0_t;
T1_t = T0_t;

%% 1--->2: compressore
P2_t = beta_c.*P1_t;
T2_ad = T1_t.*beta_c.^((gamma_a-1)/gamma_a);
T2_t = T1_t + (T2_ad-T1_t)./eta_c;

%% 2--->3: burner
P3_t = pi_b.*P2_t;
ff=(cp_a.*T2_t-cp_gc.*T3_t)./(cp_gc.*T3_t-eta_b.*deltaH);
f = [f,ff];


%% 3--->4: turbina
T4_t=T3_t-(cp_a.*(T2_t-T1_t))./(eta_c_m.*eta_t_m.*(1+ff).*cp_gc);
T4_ad=T3_t+(T4_t-T3_t)./eta_t;
P4_t=P3_t.*(T3_t./T4_ad).^(gamma_gc/(1-gamma_gc));
T5_t=T4_t;

%% hp OE
P5_st=P0;
T5_st_i=T4_t.*(P4_t./P5_st).^((1-gamma_gc)/gamma_gc);
T5_st=T4_t-eta_n.*(T4_t-T5_st_i);
a_e=sqrt(gamma_gc.*R_gc.*T5_st);
v_ee=sqrt(2*cp_gc.*(T5_t-T5_st));
M_ee=v_ee./a_e;
rho_e=P5_st./(R_gc.*T5_st);
m_tot = rho_e.*A_e.*v_ee; %aria già spillata 
m_a_a = (m_tot./(1+ff))/0.98; %considero il 2% di spillamento di portata dal compressore
m_f = ff.*m_a_a;
TH= m_tot.*v_ee - m_a_a.*v_vol;

if (M_ee > 1)
    M_ee = 1;
    T5_st = T4_t.*(2./(gamma_gc+1));
    %T5_st = (2*cp_gc*T4_t)/(gamma_gc*R_gc+2*cp_gc);
    T5_st_i = T4_t + (T5_st-T4_t)./eta_n;
    P5_st = P4_t.*(T5_st_i./T4_t).^(gamma_gc/(gamma_gc-1));
    v_ee = sqrt(2.*cp_gc.*(T4_t-T5_st));
    %controllo
    a_e = sqrt(gamma_gc.*R_gc.*T5_st);
    M_ee = v_ee./a_e;
    rho_e=P5_st./(R_gc.*T5_st);
    m_tot = rho_e.*A_e.*v_ee;
    m_a_a = (m_tot./(1+ff))/0.98; %considero il 2% di spillamento di portata dal compressore
    m_f = ff.*m_a_a;
    m_norm_ = m_a_a./(A_e.*P4_t).*sqrt(R_gc.*T4_t);
    m_norm = [m_norm,m_norm_];
    TH= m_tot.*v_ee - m_a_a.*v_vol + A_e.*(P5_st-P0);
end

THRUST = [THRUST, TH];
v_e = [v_e,v_ee];
M_e = [M_e,M_ee];

m_a = [m_a, m_a_a];
%% calcolo quantità all'ingresso della presa (considerando flusso adiabatico isoentropico)
%fun = @(M_in) m_a - (rho0_t/(1+(gamma_a-1)/2*M_in^2)^(1/(gamma_a-1))) * A_in * M_in * sqrt(gamma_a*R_a* T0_t/(1+(gamma_a-1)/2*M_in^2));
fun = @(M_in) m_a_a - P0_t.*A_in.*M_in.*sqrt(gamma_a).*(1+(gamma_a-1)./2.*M_in.^2).^(-(gamma_a+1)/(2*gamma_a-2))./sqrt(T0_t.*R_a);
%dfun = @(M_in) rho0_t*A_in* ( (1+(gamma_a-1)/2*M_in^2)^(gamma_a/(1-gamma_a))*M_in^2*(gamma_a*R_a*T1_t/(1+(gamma_a-1)/2*M_in^2))^1/2 - ...
 %   (1+(gamma_a-1)/2*M_in^2)^(1/(1-gamma_a))*(gamma_a*R_a*T1_t/(1+(gamma_a-1)/2*M_in^2))^1/2 + ...
  %  (1+(gamma_a-1)/2*M_in^2)^(1/(1-gamma_a))*M_in*1/2*(gamma_a*R_a*T1_t/(1+(gamma_a-1)/2*M_in^2))^(-1/2)*(gamma_a*R_a*T1_t*(gamma_a-1)*M_in/(1+(gamma_a-1)/2*M_in^2)^2) );

x0 = 0.3;

M_inn = fzero(fun, x0);
M_in = [M_in,M_inn];
T1 = T0_t./(1+(gamma_a-1)/2*M_inn^2);
rho1 = rho0_t./(1+(gamma_a-1)./2.*M_inn.^2).^(1/(gamma_a-1));
P1 = rho1.*R_a.*T1;
a1 = sqrt(gamma_a.*R_a.*T1);
v1 = M_inn.*a1;
m_a2 = rho1.*A_in.*v1;
v_in = [v_in; v1];
%nmax = 1000;
%toll = 1e-3;
%[xvect,it]=newton(x0,nmax,toll,fun,dfun, 1);
%M_in = xvect(end);
%% Prestazioni
P_p = TH.*v_vol;
P_dis = 1/2.*m_tot.*(v_ee-v_vol).^2;
P_av = m_f.*deltaH;
eta_p_ = P_p./(P_p+P_dis);
eta_ter_ = (P_p+P_dis)./P_av;
eta_glob_ = eta_p_.*eta_ter_;
tsfc = m_f./TH.*3600;
I_sp_ = TH./m_a_a;

eta_p = [eta_p,eta_p_];
eta_ter = [eta_ter,eta_ter_];
eta_glob = [eta_glob,eta_glob_];
TSFC = [TSFC, tsfc];
I_sp = [I_sp,I_sp_];
end
