clear
close all
clc

%aria
gamma_a=1.4;
R_a = 287;
cp_a=1004;
%gas combusti
gamma_gc=1.33;
R_gc=286.6;
cp_gc=1155;

%quota
z = 9144;
p_amb0 = 101325;
T_amb0 = 288.15;
rho_amb0 = 1.225;
a=0.0065;
g = 9.81;
T0=T_amb0-a*z;
P0=p_amb0*(1-(a*z)/T_amb0)^(g/(R_a*a));
rho0=P0/(R_a*T0);

% dati Bmw003-A-2
v_in = 108.187155101845;      %velocità presa dal ciclo a tubo di flusso variabile a 500mph
rho1 = 0.559440041203149;     %densità presa ....
v_volo = 500.*1.609/3.6;
a_volo = sqrt(gamma_a*R_a*T0);
M_volo = v_volo/a_volo;
T0_t = T0*(1+(gamma_a-1)/2*M_volo^2);
P0_t = P0*(1+(gamma_a-1)/2*M_volo^2)^(gamma_a/(gamma_a-1));
T3_t = 1026.95;     %assumo la stessa T in camera che ho a quota 0 a punto fisso
%T3=320;
%M3=0.26;
%T3=(T3-32)*5/9+273.15;
%T3_t=T3*(1+((gamma_a-1)/2)*M3^2);

% Rendimenti
% diffusore
pi_d=0.97;
% compressore
beta_c=3.09;
eta_c_tot=0.74;
eta_c=0.78;
eta_c_m=eta_c_tot/eta_c;
% burner
eta_b=0.95;
pi_b=0.95;
deltaH=43.961*10^6;
%turbina
eta_t_t = 0.78;
eta_t_m=0.95;
eta_t = eta_t_t/eta_t_m;
%ugello
eta_n=0.97;

%% Area di ingresso
D_in = 17;
D_in = D_in*2.54/100;
A_in = pi*D_in^2/4;
A_e = 0.1020;

m_a = rho1*A_in*v_in;

%% 0--->1: diffusore
P1_t = pi_d*P0_t;
T1_t = T0_t;

%% 1--->2: compressore 
P2_t = beta_c*P1_t;
T2_ad = T1_t*beta_c^((gamma_a-1)/gamma_a);
T2_t = T1_t + (T2_ad-T1_t)/eta_c;

m_a = m_a*0.98; %considero il 2% di spillamento di portata dal compressore 

%% 2--->3: burner
P3_t = pi_b*P2_t;
f=(cp_a*T2_t-cp_gc*T3_t)/(cp_gc*T3_t-eta_b*deltaH);

%% 3--->4: turbina con alimentazione LRE (150HP)
Power_lost = 111855; %potenza turbopompe LRE
eta_LRE_m = 0.95; %assunto rendimento meccanico
T4_t=T3_t-((cp_a*(T2_t-T1_t)/eta_c_m)+Power_lost/(m_a*eta_LRE_m))/(eta_t_m*(1+f)*cp_gc);
T4_id = T3_t+(T4_t-T3_t)/eta_t;
P4_t=P3_t*(T3_t/T4_id)^(gamma_gc/(1-gamma_gc));

%% ugello 
%ipotesi OE
T5_t=T4_t;
P_e = P0;
T_e_id = T4_t*(P_e/P4_t)^((gamma_gc-1)/gamma_gc);
T_e = T4_t-eta_n*(T4_t-T_e_id);
a_e=sqrt(gamma_gc*R_gc*T_e);
v_e=sqrt(2*cp_gc*(T4_t-T_e));
M_e=v_e/a_e;

%sottoespansione
if (M_e>1)
    M_e = 1; %blocco della portata
    T_e = T4_t*(2/(gamma_gc+1));
    T_e_id = T4_t - (T4_t-T_e)/eta_n;
    P_e = P4_t*(T_e_id/T4_t)^((gamma_gc)/(gamma_gc-1));
    v_e = sqrt(2*cp_gc*(T4_t-T_e));
    %check
    a_e = sqrt(gamma_gc*R_gc*T_e);
    M_e = v_e/a_e;
end

%% PRESTAZIONI

m_f = f*m_a;
m_e = m_a+m_f;
T = m_e*v_e - m_a*v_volo + (P_e-P0)*A_e %sottoespansione

I_sp = T/m_a;
TSFC = m_f/T*3600;

P_av = m_f*deltaH;
P_p = T*v_volo;
P_dis = 0.5*(m_a+m_f)*(v_e-v_volo)^2;
P_j = P_p+P_dis;
eta_th = P_j/P_av;
eta_p = P_p/P_j;
eta_g = eta_p*eta_th;

%% incrementi

T_x = 3420.44090698131;
incremento_T = (T-T_x)/T_x*100

I_sp_x = 385.921994975812;
incremento_Isp = (I_sp-I_sp_x)/I_sp_x*100

TSFC_x = 0.185614925465724;
incremento_TSFC = (TSFC-TSFC_x)/TSFC_x*100

