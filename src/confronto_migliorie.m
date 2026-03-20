clear
close all
clc

%% differenza percentuale di spinta

% spinte prese dai cicli risolutivi a velocità fissata 500mph=800km/h
T_BMW = 3189.56894944682;
T_Jumo = 3884.48912803725;

diff_T = (T_Jumo-T_BMW)/T_BMW*100

%% differenza prestazioni

eta_p_BMW = 0.613310442800823;
eta_ter_BMW = 0.152963289073819;
eta_glob_BMW = 0.0938139825541344;
TSFC_BMW = 0.195070192540889;
I_sp_BMW = 367.215930839889;

eta_p_Jumo = 0.630365435548920;
eta_ter_Jumo = 0.147950997241585;
eta_glob_Jumo = 0.0932631948160889;
TSFC_Jumo = 0.206031462052887;
I_sp_Jumo = 443.870584067178;

diff_eta_p = (eta_p_Jumo-eta_p_BMW)/eta_p_BMW*100
diff_eta_ter = (eta_ter_Jumo-eta_ter_BMW)/eta_ter_BMW*100
diff_eta_glob = (eta_glob_Jumo-eta_glob_BMW)/eta_glob_BMW*100
diff_TSFC = (TSFC_Jumo-TSFC_BMW)/TSFC_BMW*100
diff_I_sp = (I_sp_Jumo-I_sp_BMW)/I_sp_BMW*100

