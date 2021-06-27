%% do PWB at one loss
%%
EPS0=8.854187817e-12;
MU0=4e-7*pi;
c = 1/sqrt(EPS0*MU0); ETA0 = sqrt(MU0/EPS0);
%cm = 0.4; %% 1 cm = 0.4 in
cm = 0.3937; %% 1 cm = 0.4 in
%% general 
npsim=15000;
freq=linspace(75E9,110E9,npsim);

%freq = 100e9; 
sigma_wall_mu0_ori = 6e7*MU0; % all wall is the same material
mu_r = 1; 
rt=4;
app_1_S = 0.1181*0.0394*1e-2*1e-2*cm*cm; % surface area of apperture 1
app_2_S = pi*0.25*0.25*1e-2*1e-2*cm*cm; % surface area of apperture 2
app_2_S = app_2_S*rt;

pedal_S = (0.5*1.5*0.7*1.5)*1e-2*1e-2*cm*cm;

cav_1_S_solid = (2*1.5*1.5 + 4*3.5*1.5)*1e-2*1e-2*cm*cm; % in m2

cav_1_S_solid = cav_1_S_solid;
cav_2_S_solid = cav_1_S_solid; % in m2

cav_1_S = cav_1_S_solid - app_1_S - app_2_S + 2*pedal_S;
cav_2_S = cav_1_S;

for ii=1:npsim
    clear i
        a=9.1;
            invsigma_wall_mu0t(1)=(a/(2*cav_1_S_solid/(3*10^8/freq(ii))^2)*3*c/4/pi)^2/freq(ii)/mu_r*pi;

    n_loss_multipliers = length(invsigma_wall_mu0t);



    PandS_container=zeros(n_loss_multipliers,7);
    eigval_container=zeros(n_loss_multipliers,2);
    eigvec_container=zeros(n_loss_multipliers,2,2);

    ilm = 1;
    ct=1;
    for lm = invsigma_wall_mu0t

        invsigma_wall_mu0 = lm;

        sig_wall = 4*pi*cav_1_S_solid/(3*c)*sqrt(freq(ii)*mu_r*invsigma_wall_mu0/(pi)); 
        afpwb = 4*pi/(3*c)*sqrt(freq(ii)*mu_r*invsigma_wall_mu0/(pi)); 

        sigran(ct)=sig_wall;
        afpwb_rag(ct)=afpwb;
        ct=ct+1;

        %  coupling crossection cav 1 & 3 
        sig_wall_1_3 = 4*pi*cav_1_S/(3*c)*sqrt(freq(ii)*mu_r*invsigma_wall_mu0/(pi)); 
        % coupling crossection cav 2
        sig_wall_2 = 4*pi*cav_2_S/(3*c)*sqrt(freq(ii)*mu_r*invsigma_wall_mu0/(pi));

        sig_app_1_4 = app_1_S/2; % apperture 1 and 4 is the same
        sig_app_2_3 = app_2_S/2;  % apperture 2 and 3 is the same

        % make power balanced matrix 
        M = zeros(2);
        M(1,1) = sig_wall_1_3 + sig_app_1_4 + sig_app_2_3;
        M(1,2) = -sig_app_2_3;
        M(2,1) = M(1,2);
        M(2,2) = sig_wall_2 + sig_app_2_3 + sig_app_1_4;

        %
        Pinc = zeros(2,1); Pinc(1,1) = 1; 

        % solution 
        S = M\Pinc; 
        [e_val,e_vec] = eig(M);
        e_val = sort(diag(e_val)); 

        %
        PandS_container(ilm,1) = S(1); % S_1
        PandS_container(ilm,2) = S(2); % S_2
        PandS_container(ilm,3) = sig_wall_1_3*S(1); % P_W1
        PandS_container(ilm,4) = sig_wall_2*S(2); % P_W2
        PandS_container(ilm,5) = sig_app_1_4*S(1); % P_1^-
        PandS_container(ilm,6) = sig_app_2_3*(S(1)-S(2)); % P_12
        PandS_container(ilm,7) = sig_app_1_4*S(2); % P_4^-


        eigval_container(ilm,:) = e_val;

        % 
        ilm = ilm+1;
    end

    Pset(ii)=PandS_container(:,6);
end
%%
hold on
plot(freq(1:15000)/10^9,Pset);
box on
xlabel('Frequency (GHz)');
ylabel('P_{1->2} (W)');
hold off
%% only for scaled 2-cav
nom=140;   %nom here is the actual number of modes one wants to include
clear Z_rr Z_rt Z_tt Z_tr;
ldname=['F:\Shukai\RCM z\znorm_' num2str(nom) 'a1p_9d0.mat']
load(ldname);
%%
[n,m]=size(Z_rr);
clear Gr Br P1 P2 Ptst U1 U2;
for i=1:npsim
    i
    U1=rand(nom,1);
    Yra=1i*imag(Ydiag(1:nom,1:nom,i))+real(Ydiag(1:nom,1:nom,i));
    Ptst=0.5*real(U1'*Yra*U1);
    rt=Pset(ii)/Ptst;
    U1=U1*sqrt(rt);
    
    Xiaa1=reshape(Z_tt(i,:,:),nom,nom);
    Xiab1=reshape(Z_tr(i,:),nom,1);
    Xiba1=reshape(Z_rt(i,:),1,nom);
    Xibb1=reshape(Z_rr(i,:),1,1);

    Yp=1/Zavg(i);
    Yl=1/50;
     Gr=real(Ydiag(1:nom,1:nom,i));
     Br=imag(Ydiag(1:nom,1:nom,i));

    
    Ggh=Gr.^0.5;
    Yaa1=Ggh*Xiaa1*Ggh+1i.*Br;
    Ybb1=1i*imag(Yp)+real(Yp)^0.5*Xibb1*real(Yp)^0.5;   
    Yba1=Xiba1*Ggh.*sqrt(real(Yp));
    Yab1=Ggh*Xiab1.*sqrt(real(Yp));    
      
    Yeq=Yaa1-Yab1/(Ybb1+Yl)*Yba1; %% Y_L cascade
    U2(i)=-inv(Ybb1+Yl)*Yba1*U1;
   
    P1(i)=0.5*real(U1'*Yeq*U1);
    P2(i)=0.5*real(U2(i)'*Yl*U2(i));
  
end
%clear Z_rr Z_rt Z_tt Z_tr Xibb1 Xiaa1 Xiba1 Xiab1;

P29d0=P2;
U29d0=U2;
P19d0=P1;
%%
imz=abs(U2);
h=figure;
hold on;
x1=min(min(imz));
x2=max(max(imz));

dx=0.06;
x=[x1:dx:x2];
y=hist(imz,x);
yall=sum(y)*dx;
plot(x,y/yall,'.');
title('induced voltage at load, 3cav, scaled');
hold off;
imxv2_hy=x;
imyv2_hy=y/yall;
%%
h=figure;
hold on;
fig1=plot(imxv2,imyv2,'o','DisplayName','exp');
fig2=plot(imxv2_hy,imyv2_hy,'o','DisplayName','hybrid');
fig3=plot(imxv2rf_rcmf,imyv2rf_rcmf,'o','DisplayName','rcm');
%title('PWB-RCM hybrid model, 2cav induced power, scaled');
title('PWB-RCM hybrid model, 2cav induced voltage, scaled');
xlabel('Induced Voltage (V)');
%xlabel('Induced Power (W)');
ylabel('PDF');
box on;
xlim([0 1.4]);
legend;
hold off