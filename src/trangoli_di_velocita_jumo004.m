clear
close all
clc
%% JUMO004
%****DESCRIZIONE DEL PROCENDIMENTO
%I PRIMI DUE CICLI FOR ANNIDATI SERVONO PER TROVARE L'ANGOLO ALPHA DI INGRESSO AL PRIMO STADIO DEL COMPRESSORE. 
% SONO FATTI SUL RAGGIO MEDIO SU CUI DEFINIAMO UN GRADO DI REAZIONE DI 0.97
% LE MATRICI HANNO SULLE COLONNE GLI STADI E SULLE ASCISSE TUTTE LE POSSIBILI ALPHA
%IL CICLO TROVA L'ALPHA DI INGRESSO AL RAGGIO MEDIO TALE CHE IL RAPPORTO ROTORE SIA IL PIU' VICINO
%POSSIBILE AL BETA STADIO/PERDITE DI PRESSIONE


%P2TOT sono le pressioni ideali all'uscita del rotore
%P3TOT sono le pressioni reali all'uscita del rotore

ALPHA=-pi/5:pi/10000:pi/8;
%diametri e raggi
d_c_ext=21.48*0.0254;
d_c_int=d_c_ext*0.5;
n_stadi=8;

alpha1=zeros(n_stadi,length(ALPHA));
alpha2=zeros(n_stadi,length(ALPHA));
U_c=zeros(n_stadi,1);
beta1=zeros(n_stadi,length(ALPHA));
beta2=zeros(n_stadi,length(ALPHA));

om=8700*2*pi/60;


BETA_C=3.5;
BETA_st=3.5^(1/n_stadi);
R=0.9;


cp_a=1004;
gamma_a=1.4;
T1_tot=zeros(n_stadi,length(ALPHA));
P1_tot=zeros(n_stadi,length(ALPHA));
P2_tot=zeros(n_stadi,length(ALPHA));
P3_tot=zeros(n_stadi,length(ALPHA));
P1_tot_reale=zeros(n_stadi,length(ALPHA));
T3_tot=zeros(n_stadi,length(ALPHA));
T3_tot_id=zeros(n_stadi,length(ALPHA));
C1_z=zeros(n_stadi,length(ALPHA));
C2_z=zeros(n_stadi,length(ALPHA));
W1_z=zeros(n_stadi,length(ALPHA));
W2_z=zeros(n_stadi,length(ALPHA));
C1_th=zeros(n_stadi,length(ALPHA));
C2_th=zeros(n_stadi,length(ALPHA));
W1_th=zeros(n_stadi,length(ALPHA));
W2_th=zeros(n_stadi,length(ALPHA));
rapporto_rotore=zeros(n_stadi,length(ALPHA));
perdita = 0.98;%

    r_c=d_c_int/2:0.001:d_c_ext/2;
    r_m=(d_c_ext+d_c_int)/4;
    l_paletta_i=(r_c(end)-r_c(1))*100; %LUNGHEZZA DELLA PALETTA AL PRIMO STADIO DI CIRCA 10CM
    hubtotip_i=d_c_int/d_c_ext;
    r_c=[r_c,d_c_ext/2];
for j=1:length(ALPHA)
    alpha1(1,j)=ALPHA(j);
    C1=155/cos(alpha1(1,j));

    T1_tot(1,j)=2.535720558013630e+02;
    P1_tot(1,j)=4.185434568613693e+04;
    P1_tot_reale(1,j)=4.185434568613693e+04;
  
    T3_tot_vera=3.934816560120708e+02; %N.B.: data la velocità e la quota a cui scegliamo di andare, vedere ciclo camo
    eta_c=0.78;

    % FUNZIONE PER TROVARE L'ETA STADIO  
    fun=@(x) eta_c-(BETA_C^((gamma_a-1)/gamma_a)-1)/(BETA_C^((gamma_a-1)/(gamma_a*x))-1);
    eta_st0=0.8;
    eta_st=fzero(fun,eta_st0);
    f=@(y) T3_tot_vera/T1_tot(1,j)-BETA_C^((gamma_a-1)/(gamma_a*y));
    eta_st2=fzero(f,eta_st0);
   
    for i=1:n_stadi
       
        U_c(i)=r_c(end)*om;
        C1_z(i,j)=C1*cos(alpha1(i,j));
        C1_th(i,j)=C1*sin(alpha1(i,j));
        C1_vect=[C1_z(i,j);C1_th(i,j)];
        W1_vect=C1_vect-[0;U_c(i)];
        W1_z(i,j)=W1_vect(1);
        W1_th(i,j)=W1_vect(2);
        beta1(i,j)=atan(W1_th(i,j)/W1_z(i,j));

        beta2(i,j)=atan(((-R+0.5)*2*norm([0;U_c(i)]))/C1_z(i,j)-tan(alpha1(i,j)));%uso grado di reazione
        C2_z(i,j)=C1_z(i,j);%ipotesi vel assiale costante
        W2_z(i,j)=W1_z(i,j);
        W2_th(i,j)=W2_z(i,j)*tan(beta2(i,j));
        W2_vect=[W2_z(i,j);W2_th(i,j)];
        C2_vect=W2_vect+[0;U_c(i)];
        C2_th(i,j)=C2_vect(2);
        alpha2(i,j)=atan(C2_th(i,j)/C2_z(i,j));
        alpha1(i+1,j)=alpha1(i,j);%ipotesi stadio ripetuto
        

        rapporto_rotore(i,j)=(1+(om*r_c(end)*(C2_z(i,j)*tan(alpha2(i,j))-C1_z(i,j)*tan(alpha1(i,j))))/(cp_a*T1_tot(i,j)))^(gamma_a/(gamma_a-1));%eulero
        P2_tot(i,j)=rapporto_rotore(i,j)*P1_tot(i,j);
        P3_tot(i,j)=(BETA_st/perdita)*P1_tot_reale(i,j);%P3_tot sarebbe P2_tot nel caso reale
        T3_tot_id(i,j)=T1_tot(i,j)*(BETA_st^((gamma_a-1)/gamma_a));
        T3_tot(i,j)=T1_tot(i,j)+((T3_tot_id(i,j)-T1_tot(i,j))/eta_st);
        T1_tot(i+1,j)=T3_tot(i,j);
        P1_tot(i+1,j)=P2_tot(i,j);
        P1_tot_reale(i+1,j)=P3_tot(i,j)*perdita; %P1_tot_reale sarebbe P3_tot nel caso reale
       

    end
end

% IL CRITERIO PER SCEGLIERE ALPHA SELEZIONA L'ALPHA IN INGRESSO 
% ASSOCIATO AL RAPPORTO ROTORE DELL'ULTIMO STADIO PIU' VICINO
% AL BETA STADIO/PERDITE.


T1_tot=T1_tot(1:end-1,:);
alpha1=alpha1(1:end-1,:);
errore=rapporto_rotore-(BETA_st/perdita)*ones(n_stadi,length(ALPHA));
errore=errore(end,:);
[err,index]=min(errore(errore>0))%occhio a questa funzione

% QUI ENTRA IN GIOCO L'IPOTESI DI FREE VORTEX, OVVERO R*CTHETA COSTANTE
% LUNGO LA PALETTA
K1=C1_th(1,index)*r_c(end);
K2=C2_th(1,index)*r_c(end);
k=(K1+K2)/(2*om);

%_R INDICA CHE QUESTE SONO LE MATRICI VALIDE PER UN RAGGIO QUALUNQUE
%IL NUMERO DI RIGHE SONO I POSSIBILI RAGGI, MENTRE IL NUMERO DI COLONNE GLI
%STADI, ESSENDO STADIO RIPETUTO UNO STADIO VALE L'ALTRO.
alpha1_r=zeros(length(r_c),n_stadi);
alpha2_r=zeros(length(r_c),n_stadi);
U_r=zeros(length(r_c),1);
beta1_r=zeros(length(r_c),n_stadi);
beta2_r=zeros(length(r_c),n_stadi);
c1_z_r=155*ones(length(r_c),n_stadi);
c2_z_r=zeros(length(r_c),n_stadi);
w1_z_r=zeros(length(r_c),n_stadi);
w2_z_r=zeros(length(r_c),n_stadi);
c1_th_r=zeros(length(r_c),n_stadi);
c2_th_r=zeros(length(r_c),n_stadi);
w1_th_r=zeros(length(r_c),n_stadi);
w2_th_r=zeros(length(r_c),n_stadi);
c1z=155;
c2z=155;
grado_reazione_r=zeros(length(r_c),n_stadi); %MATRICE DEI GRADI DI REAZIONE
rapporto_rotore_r=zeros(length(r_c),n_stadi); %MATRICE DEI RAPPORTI ROTORE

for j=1:n_stadi

for i=1:length(r_c)
    r=r_c(i);
    U_r(i)=om*r;
    c2_th_r(i,j)=K2/r;
    c2_vect=[c2z;c2_th_r(i,j)];
    w2_vect=c2_vect-[0;U_r(i)];
    w2_z_r(i,j)=w2_vect(1);
    w2_th_r(i,j)=w2_vect(2);
    beta2_r(i,j)=atan(w2_th_r(i,j)/w2_z_r(i,j));
    alpha2_r(i,j)=atan(c2_th_r(i,j)/c2z);

    c1_th_r(i,j)=K1/r;
    alpha1_r(i,j)=atan(c1_th_r(i,j)/c1z);
    c1_vect=[c1z;c1_th_r(i,j)];
    w1_vect=c1_vect-[0;U_r(i)];
    w1_z_r(i,j)=w1_vect(1);
    w1_th_r(i,j)=w1_vect(2);
    beta1_r(i,j)=atan(w1_th_r(i,j)/w1_z_r(i,j));

    grado_reazione_r(i,j)=0.5-c1z/U_r(i)*(tan(alpha1_r(i,j))+tan(beta2_r(i,j)))/2;
    rapporto_rotore_r(i,j)=(1+(om*r*(c2z*tan(alpha2_r(i,j))-c1z*tan(alpha1_r(i,j))))/(cp_a*T1_tot(j,1)))^(gamma_a/(gamma_a-1));%eulero

end
end

%%
T1_tot_r=ones(length(r_c),n_stadi);
c1_r=zeros(length(r_c),n_stadi);
w1_r=zeros(length(r_c),n_stadi);

for i=1:n_stadi
    for j=1:length (r_c)
        T1_tot_r(j,i)=T1_tot(i,1);
    end
end

for i=1:n_stadi
    for j=1:length (r_c)
        c1_r(j,i)=sqrt((c1_th_r(j,i)).^2+(c1_z_r(j,i)).^2);
    end
end

for i=1:n_stadi
    for j=1:length (r_c)
        w1_r(j,i)=sqrt((w1_th_r(j,i)).^2+(w1_z_r(j,i)).^2);
    end
end


T1_s_r=T1_tot_r-c1_r.^2/(2*cp_a);
a_r=sqrt(gamma_a*286.6.*T1_s_r);
Mach_r=c1_r./a_r;
Mach_r_r=w1_r./a_r;

%%
% QUI PLOTTIAMO L'ANDAMENTO DELLE PRESSIONI TOTALI LUNGO I VARI STADI
% SI NOTA CHE PER IL CASO IDEALE IL SALTO DI PRESSIONE AVVIENE SONO NEL
% ROTORE E NON CI SONO PERDITE (RAPPORTO ROTORE)
% NEL CASO REALE INVECE IL SALTO DI PRESSIONE CHE AVVIENE NEL ROTORE è PIU'
% BASSO E CI SONO DELLE PERDITE NELLO STATORE (COME NELLE SLIDE DI PARAVAN)

tt = 1:17;
pressioni = [];
for i = 1:n_stadi
    pressioni = [pressioni, P1_tot(i,index)];
    pressioni = [pressioni, P2_tot(i,index)];
end
pressioni = [pressioni, P1_tot(end, index)];

pressioni_real = [];  % se il beta_stadio fosse costante
for i = 1:n_stadi
    pressioni_real = [pressioni_real, P1_tot_reale(i,index)];
    pressioni_real = [pressioni_real, P3_tot(i,index)];
end
pressioni_real = [pressioni_real, P1_tot_reale(end, index)];

pressioni_perdite = [];  % considero una perdita di pressione del 2%

figure
plot(tt, pressioni,'linewidth',2)
hold on
plot(tt, pressioni_real,'linewidth',2)
xlabel('Stadio del compressore')
ylabel('Pressione totale')
legend('pressioni','pressioni reali')

%%
% PLOTTO LA VARIAZIONE DEGLI ANGOLI IN FUNZIONE DI r/rtip LUNGO IL PRIMO STADIO
% SECONDO UN FILE CHE HO TROVATO BETA1 E BETA2 DOVREBBERO AUMENTARE MENTRE
% I NOSTRI DIMUNUISCONO, MA POTREBBE INTENDERE IL MODULO
% MI RITROVO INVECE CON I DISEGNI CARICATI DA PARAVAN


figure;
rateo=r_c./r_c(end);
plot(rateo,alpha1_r(:,1).*180/pi,'linewidth',2);
hold on
plot(rateo,beta1_r(:,1).*180/pi,'linewidth',2);
hold on 
plot(rateo,alpha2_r(:,1).*180/pi,'linewidth',2);
hold on
plot(rateo,beta2_r(:,1).*180/pi,'linewidth',2);
legend('alpha1','beta1','alpha2','beta2');
xlabel('r/r_t_i_p')
ylabel('gradi')

%PLOTTO IL GRADO DI REAZIONE IN FUNZIONE DI r/rtip (AUMENTA DALL'HUB AL TIP)
% DA NOTARE CHE IL GRADO DI REAZIONE E' DEFINITO PER UN FLUIDO
% INCOMPRIMIBILE ANCHE SULLE SLIDE DI PARAVAN
figure;
plot(rateo,grado_reazione_r(:,1),'linewidth',2)
legend('Andamento del grado di reazione lungo il raggio');
xlabel('r/r_t_i_p')
ylabel('R')
%PLOTTO IL RAPPORTO ROTORE LUNGO GLI STADI (DIMINUISCE)
figure;
X=1:0.01:n_stadi;
x=1:1:n_stadi;
Y=spline(x,rapporto_rotore_r(1,:),X);
plot(X,Y,'linewidth',2);
legend('andamento del rapporto rotore lungo gli stadi');
xlabel('stadi')
ylabel('rapporto rotore')