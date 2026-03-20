clear
close all
clc

%% dati Jumo004-B-4

% vettore delle spinte prese dal ciclo risolutivo ogni 50mph
thr= [2852.45438419876	2836.65994692953	2849.44867823506	2887.59992222184	2949.26915915711	3042.56701378415	3161.51172809154	3304.90595814932	3471.44716693421];
vel = 150:50:550;         
vv = linspace(150, 550, 401);
p4 = polyfit(vel, thr, 4);
p_val4 = polyval(p4, vv);
thr = thr.*0.4536.*9.81;

% area di ingresso
D_in = 20;
D_interno = 10; %assumiamo che l'albero motore abbia diametro (vedi foto)
D_in = D_in*2.54/100;
D_interno = D_interno*2.54/100;
A_in = pi*D_in^2/4 - pi*D_interno^2/4;

% diffusore
pi_d=0.97;

% compressore
beta_c=3.5;
eta_c_m=0.95; 
eta_c=0.78;

% burner
eta_b=0.95;
pi_b=0.95;
deltaH=41.868*10^6;
T3_intermedia = 754.40+273.15; %T in camera di combustione

% after-burner (ulteriore iniezione combustibile prima della turbina)
T3_t = 870+273.15;

%turbina
eta_t = 0.795;
eta_t_m=0.95;

%ugello
eta_n=0.97;

% area di efflusso
A_e = 0.1020; %si considera costante oltre 8000m (preso dal BMW003)


%% dati quota 30000ft

%aria
gamma_a=1.4;
R_a = 287;
cp_a=1004;

%gas combusti
gamma_gc=1.33;
R_gc=286.6;
cp_gc=1155;

%after-burner
gamma_AB = 1.30;
cp_AB = 1199;
R_AB = 286.6;

%quota
z = 9144;
p_amb0 = 101325;
T_amb0 = 288.15;
rho_amb0 = 1.225;
a=0.0065;
g = 9.81;
T0=T_amb0-a*z;
P0=p_amb0*(1-(a*z)/T_amb0).^(g/(R_a*a));
rho0=P0/(R_a*T0);


%% risoluzione ciclo termodinamico

% inizializzazione vettori
THRUST = [];
m_a = [];
v_in = [];
v_e = [];
M_e = [];
m_norm = [];
m_f_AB = [];
f_AB = [];
f = [];
eta_p = [];
eta_ter = [];
eta_glob = [];
TSFC_AB = [];
I_sp_AB = [];
M_in = [];
A0 = [];
v_comp = [];

% ciclo per ogni velocità 
v_volo = 150:50:550; 
for i = 1:length(v_volo)
v_vol = v_volo(i).*1.609/3.6;
a_volo = sqrt(gamma_a*R_a*T0);
M_volo = v_vol./a_volo;
T0_t = T0.*(1+(gamma_a-1)./2.*M_volo.^2);
P0_t = P0.*(1+(gamma_a-1)./2.*M_volo.^2).^(gamma_a/(gamma_a-1));
rho0_t = rho0.*(1+(gamma_a-1)./2.*M_volo.^2).^(1/(gamma_a-1));

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
f_AB_ = ((1+ff).*cp_gc.*T3_intermedia-(1+ff).*cp_AB.*T3_t)./(cp_AB.*T3_t-eta_b.*deltaH); %combustione AB

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
M_e = [M_e,M_ee];
m_a = [m_a, m_a_a];
m_f_AB = [m_f_AB,m_f_AB_];
f_AB = [f_AB,f_AB_];
f = [f,ff];

%% tubo di flusso di cattura (considerando flusso adiabatico isoentropico)
A0_ = m_a_a./(rho0.*v_vol);
A0 = [A0,A0_];
a0 = sqrt(gamma_a*R_a*T0);
M0 = v_vol./a0;

P1_tt = P0_t;
T1_tt = T0_t;

v0 = v_vol;
fun = @(v1) log(A_in/A0_) + log(v1/v0) + ...
            cp_a/(gamma_a*R_a)*log((gamma_a*R_a*T1_tt- ...
            gamma_a*R_a*v1^2/(2*cp_a))/(gamma_a*R_a*T1_tt- ...
            gamma_a*R_a*v0^2/(2*cp_a)));
v_guess = 100;
v1 = fzero(fun, v_guess);
T1 = T1_tt - v1^2/(2*cp_a);
a1 = sqrt(gamma_a*R_a*T1);
M1 = v1/a1;
P1 = P1_tt/(1+(gamma_a-1)/2*M1^2)^(gamma_a/(gamma_a-1));
rho1 = P1/(R_a*T1);

v_in = [v_in,v1];
M_in = [M_in,M1];

%% presa d'aria
Pt_presa = pi_d.*P1_tt;
Tt_presa = T1_tt;
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

eta_p_x = [0.230310063207635	0.291603130591917	0.345918380121536	0.397438434989330	0.458417245250604	0.518021917729429	0.575787977787511	0.631167178425285	0.683605619373668];
incremento_eta_p = (eta_p-eta_p_x)./eta_p_x*100

eta_ter_x = [0.145557173440220	0.149054845000647	0.153457323113483	0.157338140789678	0.157013769626423	0.157231897547504	0.158024542633702	0.159395289694612	0.161322297443996];
incremento_eta_ter = (eta_ter-eta_ter_x)./eta_ter_x*100

eta_glob_x = [0.0335232818153418	0.0434648594320817	0.0530837086292029	0.0625322244395803	0.0719778197385578	0.0814495690957948	0.0909886318438556	0.100605075250829	0.110280829062986];
incremento_eta_glob = (eta_glob-eta_glob_x)./eta_glob_x*100

%% grafici
vel2 = vel.*1.609./3.6;
vv_ms = vv.*1.609./3.6;
p_val_N = p_val4;
p4_th = polyfit(vel2,THRUST,4);
p_val_th = polyval(p4_th,vv_ms);
figure
plot(vv_ms,p_val_N,'LineWidth',2)
xlabel('Velocità di volo [m/s]')
ylabel('Spinta [N]')
hold on 
grid on
plot(vv_ms,p_val_th,'LineWidth',2)
legend('AB spento','AB acceso')
