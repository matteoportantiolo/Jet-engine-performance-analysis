clear
close all
clc

% dati Jumo004
z = 9144;    %quota
thr= [680; 675; 680; 690; 700; 720; 740; 760; 790];    % vettore delle spinte prese dal grafico ogni 50mph
vel = 150:50:550;
v_volo = 150:50:550;          % subsonico in efflusso fino a 320; la portata la calcolo in efflusso
vv = linspace(150, 550, 401);
p4 = polyfit(vel, thr, 4);
p_val4 = polyval(p4, vv);

thr = thr.*0.4536.*9.81;
thr = thr';
% Rendimenti
% diffusore
pi_d=0.97; %assunto

% compressore
beta_c=3.5; 
eta_c_m=0.95; 
eta_c=0.78;

% burner
eta_b=0.95;
pi_b=0.95; %assunto
deltaH=41.868*10^6;
%turbina
eta_t_m=0.95; %assunto
eta_t = 0.795; 
%ugello
eta_n=0.97; %assunto

[THRUST, TSFC, T0, P0, rho0, T0_t, P0_t, rho0_t, A_in, A0, T1, P1, rho1, a1, v_in, v_comp, M_volo, M_in, M_e, ...
                v_e, f, m_a, m_f, m_tot, P_p, P_dis, P_av, eta_p, eta_ter, eta_glob, P5_st, m_norm,I_sp,rho_e,T5_st,P1_tt,T1_tt,P2_tt,T2_tt] = Jumo004 (z, ...
                v_volo, pi_d, beta_c, eta_c_m, eta_c, eta_b, pi_b, deltaH, eta_t, eta_t_m, eta_n);
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

