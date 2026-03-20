clear
close all
clc

%% dati

d_in = 0.508;
Q = 41.868e6;

%pressioni
beta_c = 3.5;
pi_b = 0.95; %assunto
pi_d = 0.97; %assunto

%rendimenti
eta_c = 0.78;
eta_b = 0.95;
eta_t = 0.795;
eta_m = 0.97; %assunto
eta_n = 0.97; %assunto
T_4 = 754.40+273.15;

%aria
gamma_a = 1.4;
R_a = 287;
cp_a = 1004;

%gas combusti
gamma_gc = 1.33;
R_gc = 286.58;
cp_gc = 1155;

%quota
z = 9144;
p_amb0 = 101325;
T_amb0 = 288.15;
rho_amb0 = 1.225;
a=0.0065;
g = 9.81;
T_amb=T_amb0-a*z;
p_amb=p_amb0*(1-(a*z)/T_amb0)^(g/(R_a*a));
rho_amb=p_amb/(R_a*T_amb);

%% risoluzione

% INGRESSO

v_in = 153.5980;
v_in = v_in*0.44704;
a_in = sqrt(gamma_a*R_a*T_amb);
M_in = v_in/a_in;
m_a = rho_amb*((d_in^2*pi)/4)*v_in;

T_1tot = T_amb*(1+((gamma_a-1)/2)*M_in^2);
p_1tot = p_amb*(1+((gamma_a-1)/2)*M_in^2)^(gamma_a/(gamma_a-1));

% DIFFUSORE

p_2tot = pi_d*p_1tot;
T_2tot = T_1tot;

% COMPRESSORE

p_3tot = beta_c*p_2tot;
T_3tot_id = T_2tot*((beta_c)^((gamma_a-1)/gamma_a));
T_3tot = T_2tot + ((T_3tot_id-T_2tot)/eta_c);

% BURNER

p_4tot = pi_b*p_3tot;

M_b = 0.15; %assunto dal BMW003
T_4tot = T_4*(1+((gamma_a-1)/2)*M_b^2);
f = ((cp_gc*T_4tot)-(cp_a*T_3tot))/((Q*eta_b)-(cp_gc*T_4tot));
m_f = f*m_a;

% TURBINA

T_5tot = T_4tot - (m_a*cp_a*(T_3tot-T_2tot)/eta_m)/(eta_m*(m_a+m_f)*cp_gc);
T_5tot_id = T_4tot + ((T_5tot-T_4tot)/eta_t);
p_5tot = p_4tot*((T_5tot_id/T_4tot)^((gamma_gc)/(gamma_gc-1)));

% AFTER-BURNER

%dati
Q_AB = Q;
f_AB = f; %T = 870+273.15
m_fAB = f_AB*m_a;
eta_AB = eta_b;
pi_AB = pi_b;
gamma_AB = 1.30;
cp_AB = 1243;
R_AB = 286.6;

T_6tot = ((1+f)*cp_gc*T_5tot+f_AB*Q_AB*eta_AB)/((1+f+f_AB)*cp_AB);
p_6tot = pi_AB*p_5tot;

if p_6tot<p_amb
    fprintf('errore: espansione in ugello non possibile');
end

% UGELLO

T_7tot = T_6tot;

%ipotesi OE
p_7 = p_amb;
T_7_id = T_6tot*(p_7/p_6tot)^((gamma_AB-1)/(gamma_AB));
T_7 = T_6tot - eta_n*(T_6tot-T_7_id);
a_e = sqrt(gamma_AB*R_AB*T_7);
v_e = sqrt(2*cp_AB*(T_6tot-T_7));
M_e = v_e/a_e; %ipotesi giusta

%sottoespansione
if (M_e>1)
    M_e = 1; %blocco della portata
    T_7 = T_6tot*(2/(gamma_AB+1));
    T_7_id = T_6tot - (T_6tot-T_7)/eta_n;
    p_7 = p_6tot*(T_7_id/T_6tot)^((gamma_AB)/(gamma_AB-1));
    v_e = sqrt(2*cp_AB*(T_6tot-T_7));
    %check
    a_e = sqrt(gamma_AB*R_AB*T_7);
    M_e = v_e/a_e;
end

%check2
T_7tot_check = T_7*(1+(gamma_AB-1)/2);
p_7tot_check = p_7*(1+(gamma_AB-1)/2)^(gamma_AB/(gamma_AB-1));

% SPINTA

m_e = m_a+m_f+m_fAB;
p_e = p_7;
T_e = T_7;
rho_e = p_e/(R_AB*T_e);
A_e = m_e/(rho_e*v_e);

T_AB = m_e*v_e - m_a*v_in + (p_e-p_amb)*A_e; 
T_AB/4.44822  %T_x = 680;

T_x = 680; %variabile
T_x = T_x*4.44822;
incremento_T = (T_AB-T_x)/T_x*100

% PRESTAZIONI

I_sp_AB = T_AB/m_a
I_sp = 4.286655975858299e+02;
incremento_Isp = (I_sp_AB-I_sp)/I_sp*100

TSFC_AB = (m_f+m_fAB)/T_AB*3600
TSFC = 0.180959084752399;
incremento_TSFC = (TSFC_AB-TSFC)/TSFC*100

P_av = m_f*Q+m_fAB*Q_AB;
P_p = T_AB*v_in;
P_dis = 0.5*m_e*(v_e-v_in)^2;
P_j = P_p+P_dis;
eta_th = P_j/P_av;
eta_p = P_p/P_j;
eta_g = eta_p*eta_th;

