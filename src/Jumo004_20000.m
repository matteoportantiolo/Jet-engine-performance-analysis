clear
close all
clc

%% dati

d_in = 0.508;
T_x = 4448.22;
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
z = 6096;
p_amb0 = 101325;
T_amb0 = 288.15;
rho_amb0 = 1.225;
a=0.0065;
g = 9.81;
T_amb=T_amb0-a*z;
p_amb=p_amb0*(1-(a*z)/T_amb0)^(g/(R_a*a));
rho_amb=p_amb/(R_a*T_amb);

%% risoluzione

v_in = [70:0.01:80];
i = 0;
T = 0;

while (abs(T-T_x)>10)
i = i+1;

% INGRESSO

a_in = sqrt(gamma_a*R_a*T_amb);
M_in = v_in(i)/a_in;
m_a = rho_amb*((d_in^2*pi)/4)*v_in(i);

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

% UGELLO

T_6tot = T_5tot;

%ipotesi OE
p_6 = p_amb;
T_6_id = T_5tot*(p_6/p_5tot)^((gamma_gc-1)/(gamma_gc));
T_6 = T_5tot - eta_n*(T_5tot-T_6_id);
a_e = sqrt(gamma_gc*R_gc*T_6);
v_e = sqrt(2*cp_gc*(T_5tot-T_6));
M_e = v_e/a_e; %ipotesi sbagliata

%sottoespansione
if (M_e>1)
    M_e = 1; %blocco della portata
    T_6 = T_5tot*(2/(gamma_gc+1));
    T_6_id = T_5tot - (T_5tot-T_6)/eta_n;
    p_6 = p_5tot*(T_6_id/T_5tot)^((gamma_gc)/(gamma_gc-1));
    v_e = sqrt(2*cp_gc*(T_5tot-T_6));
    %check
    a_e = sqrt(gamma_gc*R_gc*T_6);
    M_e = v_e/a_e;
end

%check2
T_6tot_check = T_6*(1+(gamma_gc-1)/2);
p_6tot_check = p_6*(1+(gamma_gc-1)/2)^(gamma_gc/(gamma_gc-1));

% PRESTAZIONI

m_e = m_a+m_f;
p_e = p_6;
T_e = T_6;
rho_e = p_e/(R_gc*T_e);
A_e = m_e/(rho_e*v_e);

T = m_e*v_e - m_a*v_in(i) + (p_e-p_amb)*A_e;
TSFC = m_f/T*3600;
end

v_in(i)*3.6 %km/h
v_in(i)*2.2369 %mph
