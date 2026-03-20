clear
close all
clc

% dati Bmw003-A-2
z = 9144;    %quota
thr= [635; 620; 620; 630; 650; 670; 710; 750; 800];    % vettore delle spinte prese dal grafico ogni 50mph
vel = 150:50:550;
v_volo = 500;          % velocità a cui accendere LRE
vv = linspace(150, 550, 401);
p4 = polyfit(vel, thr, 4);
p_val4 = polyval(p4, vv);
flag = 0;
thrust = zeros(length(v_volo),1);

for j = 1:length(v_volo)
    i = 1;  
    while i <= length(vv) & flag==0
        if(vv(i) == v_volo(j))
            thrust(j) = p_val4(i);         % trovo la thrust alla velocità di volo che ho scelto dal polinomio approssimante
            flag = 1;
        end    
        i = i+1;
    end
end

thrust=thrust.*0.4536.*9.81;
thr = thr.*0.4536.*9.81;
thr = thr';
% Rendimenti
% diffusore
pi_d=0.97;

% compressore
beta_c=3.09;
eta_c_tot=0.74;
eta_c=0.78;

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

[THRUST, TSFC, T0, P0, rho0, T0_t, P0_t, rho0_t, A_in, A0, T1, P1, rho1, a1, v_in, v_comp, M_volo, M_in, M_e, ...
                v_e, f, m_a, m_f, m_tot, P_p, P_dis, P_av, eta_p, eta_ter, eta_glob, P5_st,I_sp] = BMW003_LRE (z, ...
                v_volo, pi_d, beta_c, eta_c_tot, eta_c, eta_b, pi_b, deltaH, eta_t_t, eta_t_m, eta_n);
THRUST
thr;
