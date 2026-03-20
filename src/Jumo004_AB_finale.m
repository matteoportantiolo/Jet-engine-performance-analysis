clear
close all
clc

% dati Jumo004
z = 9144;    %quota
thr= [2852.45438419876	2836.65994692953	2849.44867823506	2887.59992222184	2949.26915915711	3042.56701378415	3161.51172809154	3304.90595814932	3471.44716693421];    % vettore delle spinte prese dal ciclo
vel = 150:50:550;
v_volo = 150:50:550;          % subsonico in efflusso fino a 320; la portata la calcolo in efflusso
vv = linspace(150, 550, 401);
p4 = polyfit(vel, thr, 4);
p_val4 = polyval(p4, vv);


thr = thr';
% Rendimenti
% diffusore
pi_d=0.97; %assunto

% compressore
beta_c=3.5; 
eta_c_m=0.95; %assunto
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

[THRUST, TSFC_AB, T0, P0, rho0, T0_t, P0_t, rho0_t, A_in, A0, T1, P1, rho1, a1, v_in, v_comp, M_volo, M_in, M_e, ...
                v_e, f, f_AB,m_a, m_f,m_f_AB, m_tot, P_p, P_dis, P_av, eta_p, eta_ter, eta_glob, P5_st, m_norm,I_sp_AB,rho_e,T5_st] = Jumo004_AB (z, ...
                v_volo, pi_d, beta_c, eta_c_m, eta_c, eta_b, pi_b, deltaH, eta_t, eta_t_m, eta_n);
THRUST
thr
f_AB

%% incrementi

T = [2852.45438419876	2836.65994692953	2849.44867823506	2887.59992222184	2949.26915915711	3042.56701378415	3161.51172809154	3304.90595814932	3471.44716693421];
incremento_T = (THRUST-T)./T*100

I_sp = [433.757851130652	420.798815493707	410.046894641573	401.201760827860	394.028173467320	388.636284213958	384.172871406587	380.281354468354	376.626461366534];
incremento_Isp = (I_sp_AB-I_sp)./I_sp*100

TSFC = [0.180991710390398	0.185988751659193	0.190104515598963	0.193344951595921	0.195720909906039	0.197097960809693	0.197853466789393	0.198145656554698	0.198135073812149];
incremento_TSFC = (TSFC_AB-TSFC)./TSFC*100

%% grafici
vel2 = vel.*1.609./3.6;
vv_ms = vv.*1.609./3.6;
p_val_N = p_val4;
p4_th = polyfit(vel2,THRUST,4);
p_val_th = polyval(p4_th,vv_ms);
figure
plot(vv_ms,p_val_N)
xlabel('velocità di volo[m/s]')
ylabel('spinta[N]')
hold on 
grid on
plot(vv_ms,p_val_th)
legend('AB spento','AB acceso')


