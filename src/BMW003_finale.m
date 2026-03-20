clear
close all
clc

% dati Bmw003-A-2
z = 9144;    %quota
thr= [635; 620; 620; 630; 650; 670; 710; 750; 800];    % vettore delle spinte prese dal grafico ogni 50mph
vel = 150:50:550;
v_volo = 150:50:550;          % subsonico in efflusso fino a 320; la portata la calcolo in efflusso
vv = linspace(150, 550, 401);
p4 = polyfit(vel, thr, 4);
p_val4 = polyval(p4, vv);

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
deltaH=43.961*10^6;
%turbina
eta_t_t = 0.78;
eta_t_m=0.95;
eta_t = eta_t_t/eta_t_m;
%ugello
eta_n=0.97;

[THRUST, TSFC, T0, P0, rho0, T0_t, P0_t, rho0_t, A_in, T1, P1, rho1, a1, v_in,v_comp, M_volo, M_in, M_e, ...
                v_e, f, m_a, m_f, m_tot, P_p, P_dis, P_av, eta_p, eta_ter, eta_glob,P5_st, m_norm,I_sp,P1_tt,T1_tt,P2_tt,T2_tt] = BMW003 (z, ...
                v_volo, pi_d, beta_c, eta_c_tot, eta_c, eta_b, pi_b, deltaH, eta_t_t, eta_t_m, eta_n);
THRUST
thr

%calcolo errore
err = [];
for i=1:length(thr)
    err(i) = abs(THRUST(i)-thr(i))/thr(i)*100;
end
err

%grafici
vel2 = vel.*1.609./3.6;
vv_ms = vv.*1.609./3.6;
p_val_N = p_val4.*0.4536.*9.81;
p4_th = polyfit(vel2,THRUST,4);
p_val_th = polyval(p4_th,vv_ms);
figure
plot(vv_ms,p_val_N)
xlabel('velocità di volo[m/s]')
ylabel('spinta[N]')
hold on 
grid on
plot(vv_ms,p_val_th)
legend('Dal grafico','Dai calcoli')
