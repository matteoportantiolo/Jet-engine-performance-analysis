function [THRUST, TSFC_AB, T0, P0, rho0, T0_t, P0_t, rho0_t, A_in, A0, T1, P1, rho1, a1, v_in, v_comp, M_volo, M_in, M_e, ...
                v_e, f, f_AB,m_a, m_f,m_f_AB, m_tot, P_p, P_dis, P_av, eta_p, eta_ter, eta_glob, P5_st, m_norm,I_sp_AB,rho_e,T5_st]...
                = Jumo004_AB (z, v_volo, pi_d, beta_c, eta_c_m, eta_c, eta_b, pi_b, deltaH, eta_t, eta_t_m, eta_n)

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
m_f_AB = [];
f_AB = [];
f = [];
M_e = [];
m_norm = [];
eta_p = [];
eta_ter = [];
eta_glob = [];
TSFC_AB = [];
I_sp_AB = [];
M_in = [];
A0 = [];
v_comp = [];

for i = 1:length(v_volo)
v_vol = v_volo(i).*1.609/3.6;
a_volo = sqrt(gamma_a*R_a*T0);
M_volo = v_vol./a_volo;
T0_t = T0.*(1+(gamma_a-1)./2.*M_volo.^2);
P0_t = P0.*(1+(gamma_a-1)./2.*M_volo.^2).^(gamma_a/(gamma_a-1));
rho0_t = rho0.*(1+(gamma_a-1)./2.*M_volo.^2).^(1/(gamma_a-1));
T3_intermedia = 754.40+273.15; %T in camera di combustione
T3_t = 870+273.15; %T in camera di combustione con AB che si alza,
                   %l'AB di fatto e posto prima della turbina
                   %in serie al combustore
gamma_AB = 1.30;
cp_AB_vero = 1243;
cp_AB = (cp_AB_vero+cp_gc)/2;
R_AB = 286.6;

%eta_d = ((1+(gamma_a-1)/2*M_volo^2)*pi_d^((gamma_a-1)/gamma_a)-1)/((gamma_a-1)/2*M_volo^2);    %non l'ho usato; 
% serve se si considera il flusso come non adiabatico (a velocità basse si
% abbassa molto (guardare slide 14 della presa d'aria)

%% Area di ingresso
D_in = 20;
D_interno = 10; %assumiamo che l'albero motore abbia diametro (vedi foto)
D_in = D_in*2.54/100;
D_interno = D_interno*2.54/100;
A_in = pi*D_in^2/4 - pi*D_interno^2/4;

A_e = 0.1020;      % la considero costante tanto varia al massimo tra 970 e 1020 cm^2
                 % presa dal BMW003

%% 0--->1: diffusore
P1_t = pi_d.*P0_t;
T1_t = T0_t;

%% 1--->2: compressore
P2_t = beta_c.*P1_t;
T2_ad = T1_t.*beta_c.^((gamma_a-1)/gamma_a);
T2_t = T1_t + (T2_ad-T1_t)./eta_c;

%% 2--->3: burner
P3_t = pi_b.*P2_t;
ff=(cp_a.*T2_t-cp_gc.*T3_intermedia)./(cp_gc.*T3_intermedia-eta_b.*deltaH); %combustione burner
f_AB_ = (cp_gc.*T3_intermedia-cp_AB.*T3_t)./(cp_AB.*T3_t-eta_b.*deltaH); %combustione AB

%% 3--->4: turbina
T4_t=T3_t-(cp_a.*(T2_t-T1_t))./(eta_c_m.*eta_t_m.*(1+ff).*cp_AB);
T4_ad=T3_t+(T4_t-T3_t)./eta_t;
P4_t=P3_t.*(T3_t./T4_ad).^(gamma_AB/(1-gamma_AB));
T5_t=T4_t;

%% 4--->5: ugello
% ipotesi OE
P5_st=P0;
T5_st_i=T4_t.*(P4_t./P5_st).^((1-gamma_AB)/gamma_AB);
T5_st=T4_t-eta_n.*(T4_t-T5_st_i);
a_e=sqrt(gamma_AB.*R_AB.*T5_st);
v_ee=sqrt(2*cp_AB.*(T5_t-T5_st));
M_ee=v_ee./a_e;
rho_e=P5_st./(R_AB.*T5_st);
m_tot = rho_e.*A_e.*v_ee; %aria già spillata 
m_a_a = (m_tot./(1+ff+f_AB_))/0.98; %considero il 2% di spillamento di portata dal compressore
m_f = ff.*m_a_a; 
m_f_AB_ = f_AB_.*m_a_a;
TH= m_tot.*v_ee - m_a_a.*v_vol;

% sottoespansione
if (M_ee > 1)
    M_ee = 1;
    T5_st = T4_t.*(2./(gamma_AB+1));
    T5_st_i = T4_t + (T5_st-T4_t)./eta_n;
    P5_st = P4_t.*(T5_st_i./T4_t).^(gamma_AB/(gamma_AB-1));
    v_ee = sqrt(2.*cp_AB.*(T4_t-T5_st));
    a_e = sqrt(gamma_AB.*R_AB.*T5_st);
    M_ee = v_ee./a_e; %check
    rho_e=P5_st./(R_AB.*T5_st);
    m_tot = rho_e.*A_e.*v_ee;
    m_a_a = (m_tot./(1+ff+f_AB_))/0.98; %considero il 2% di spillamento di portata dal compressore
    m_f = ff.*m_a_a;
    m_f_AB_ = f_AB_.*m_a_a;
    m_norm_ = m_a_a./(A_e.*P4_t).*sqrt(R_AB.*T4_t);
    m_norm = [m_norm,m_norm_];
    TH= m_tot.*v_ee - m_a_a.*v_vol + A_e.*(P5_st-P0);
end

THRUST = [THRUST, TH];
v_e = [v_e,v_ee];
m_f_AB = [m_f_AB,m_f_AB_];
f_AB = [f_AB,f_AB_];
f = [f,ff];
M_e = [M_e,M_ee];

m_a = [m_a, m_a_a];
%% tubo di flusso di cattura (considerando flusso adiabatico isoentropico)
A0_ = m_a_a./(rho0.*v_vol);
A0 = [A0,A0_];
a0 = sqrt(gamma_a*R_a*T0);
M0 = v_vol./a0;
T0_t = T0.*(1+(gamma_a-1)./2.*M0.^2);
P0_t = P0.*(1+(gamma_a-1)./2.*M0.^2).^(gamma_a/(gamma_a-1));
% assumiamo flusso adiabatico isoentropico prima della presa
P1_t = P0_t;
T1_t = T0_t;

v0 = v_vol;
fun = @(v1) log(A_in/A0_) + log(v1/v0) + ...
            cp_a/(gamma_a*R_a)*log((gamma_a*R_a*T1_t-gamma_a*R_a*v1^2/(2*cp_a))/(gamma_a*R_a*T1_t-gamma_a*R_a*v0^2/(2*cp_a)));
v_guess = 100;
v1 = fzero(fun, v_guess);
T1 = T1_t - v1^2/(2*cp_a);
a1 = sqrt(gamma_a*R_a*T1);
M1 = v1/a1;
P1 = P1_t/(1+(gamma_a-1)/2*M1^2)^(gamma_a/(gamma_a-1));
rho1 = P1/(R_a*T1);

v_in = [v_in,v1];
M_in = [M_in,M1];

%% presa d'aria
Pt_presa = pi_d.*P1_t;
Tt_presa = T1_t;
M2 = 0.3; %assunto
P2 = Pt_presa./(1+(gamma_a-1)/2.*M2.^2).^(gamma_a/(gamma_a-1));
T2 = Tt_presa./(1+(gamma_a-1)/2.*M2.^2);
a2 = sqrt(gamma_a*R_a.*T2);
v2 = M2.*a2;
v_comp = [v_comp,v2];

%% Prestazioni
P_p = TH.*v_vol;
P_dis = 1/2.*m_tot.*(v_ee-v_vol).^2;
P_av = (m_f+m_f_AB_).*deltaH;
eta_p_ = P_p./(P_p+P_dis);
eta_ter_ = (P_p+P_dis)./P_av;
eta_glob_ = eta_p_.*eta_ter_;
tsfc = (m_f+m_f_AB_)./TH.*3600;
I_sp_AB_ = TH./m_a_a;

eta_p = [eta_p,eta_p_];
eta_ter = [eta_ter,eta_ter_];
eta_glob = [eta_glob,eta_glob_];
TSFC_AB = [TSFC_AB, tsfc];
I_sp_AB = [I_sp_AB,I_sp_AB_];
end
