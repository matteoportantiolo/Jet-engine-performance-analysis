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

thrust=700;
thrust=thrust*4.44822;
v_vect = [350:0.01:500];
T = [];
j = 0;
T(1) = 0;

for i = 1:length(v_vect)
j = j+1;

% dati Bmw003-A-2
v_in = v_vect(i)*0.44704;
a_in = sqrt(gamma_a*R_a*T0);
M_in = v_in/a_in;
T0_t = T0*(1+(gamma_a-1)/2*M_in^2);
P0_t = P0*(1+(gamma_a-1)/2*M_in^2)^(gamma_a/(gamma_a-1));
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
deltaH=41.868*10^6;
%turbina
eta_t_t = 0.78;
eta_t_m=0.95;
eta_t = eta_t_t/eta_t_m;
%ugello
eta_n=0.97;

%% Area di ingresso
D_in = 20;
D_in = D_in*2.54/100;
A_in = pi*D_in^2/4;

%% 0--->1: diffusore
P1_t = pi_d*P0_t;
T1_t = T0_t;

%% 1--->2: compressore
P2_t = beta_c*P1_t;
T2_ad = T1_t*beta_c^((gamma_a-1)/gamma_a);
T2_t = T1_t + (T2_ad-T1_t)/eta_c;

%% 2--->3: burner
P3_t = pi_b*P2_t;
f=(cp_a*T2_t-cp_gc*T3_t)/(cp_gc*T3_t-eta_b*deltaH);


%% 3--->4: turbina
T4_t=T3_t-(cp_a*(T2_t-T1_t))/(eta_c_m*eta_t_m*(1+f)*cp_gc);
T4_ad=T3_t+(T4_t-T3_t)/eta_t;
P4_t=P3_t*(T3_t/T4_ad)^(gamma_gc/(1-gamma_gc));
T5_t=T4_t;

%% hp OE
P5_st=44551;
T5_st_i=T4_t*(P4_t/P5_st)^((1-gamma_gc)/gamma_gc);
T5_st=T4_t-eta_n*(T4_t-T5_st_i);
a_e=sqrt(gamma_gc*R_gc*T5_st);
v_e=sqrt(2*cp_gc*(T4_t-T5_st));
M_e=v_e/a_e;
rho_e=P5_st/(R_gc*T5_st);
%A_e=(m_a+m_f)/(v_e*rho_e);
m_a = rho0*A_in*v_in;
m_f = f*m_a;
T(j) = (m_a+m_f)*v_e-m_a*v_in;
end

for j = 1:length(T)
    if abs(T(j)-thrust) < 0.05
        v_vect(j) %mph
        v_vect(j)*0.44704 %m/s
        v_vect(j)*3.6*0.44704 %km/h
    end
end



