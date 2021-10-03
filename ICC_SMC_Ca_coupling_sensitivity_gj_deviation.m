function [T,Y] = ICC_SMC_Ca_coupling_sensitivity_gj_deviation(number,Mu)
clc, close all
warning('off')

%% Define simulation time
TMAX  = 900;   %[ms] simulation time
tspan = [0 TMAX];

%% Define connectivity (index) matrix
parent = 0:number-1;
n = length(parent);
m = find(parent);
Index = sparse (parent(m), m, 1, n, n);
Index = full(Index);
Index = Index+triu(Index,1).';


%% constants
Cm = .025;
G_max_Na = 20, E_Na = 55.21;
Temp = 310, T_exp = 297;
F = 96.4846, R = 8.3144;
Q10Na = 2.45, Q10K = 1.5, Q10Ca = 2.1;
G_max_kv = 6.3, E_K = -75.91;
G_max_ERG = 2.5, G_max_bk = .15;

G_max_CaCl = 10.1, E_Cl = -11.23;
G_max_Ltyp = 2;
fc_ICC = 0.01;
Vol_ICC = 1e-12, P_cyto_ICC = 0.7, P_PU = 0.001;
P_mito_PU = 0.12871, P_ER_PU = .1;
Jmax_NaCa = 0.05, b_PU = .5;
Cmito_PU = 0.006995;
g_H_PU = 0.0033333, deltapH_PU = -0.4;
rho_res_PU = 0.4, g_PU = 0.85;

K_Na_PU = 9.4, 
Na_i = 30, 
n_PU = 2, K_Ca_PU = 0.003;

ra_PU = 6.394e-10, rb_PU = 1.762e-13;
rc1_PU = 2.656e-19, rc2_PU = 8.632e-27;
r1_PU = 2.077e-18, r2_PU = 1.728e-9, r3_PU = 1.059e-26;

deltaPsi_B = 50;
Kres_PU = 1.35e18;
tot_NAD_m = 8;
J_red_bas = 0.3333;
Jmax_leaK = 0.01;
u1_PU = 15, u2_PU = 1.1;
KCa_PDH_PU = 0.00005;

beta_max = 2.055,beta1_PU = 1.66,beta2_PU  =0.0249, beta3_PU = 4, beta4_PU = 2.83, beta5_PU = 1.3, beta6_PU = 2.66, beta7_PU = 0.16;
Glc_PU = 1;
tot_ANP_i = 2, tot_ANP_m = 12;
Jmax_ANT = 15;
rho_F1_PU = 0.7;

pa_PU = 1.656e-5, pb_PU = 3.373e-7, pc1_PU = 9.651e-14, pc2_PU = 4.845e-19;
p1_PU = 1.346e-8, p2_PU = 7.739e-7,p3_PU = 6.65e-15;
K_F1_PU = 1.71e9, Pi_m_PU = 20;
frac_PU = 0.5;

Jmax_uni = 5000, deltaPsi_s = 91;
conc_PU = 0.001, K_trans_PU = 0.006;
L_PU = 50,K_act_PU = 0.00038, Na_PU = 2.8;

Jmax_IP3 = 50000, d_ACT_PU = 0.001, J_ERleak = 1.666667, d_IP3_PU = 0.00025;
Jmax_serca = 1.8333 , K_serca_PU = 0.00042;
d_INH_PU = 0.0014, tauh_PU = 4;
fe_PU = .01, fm_PU = .0003;
J_max_PMCA = 0.088464;

Ca_o = 2.5;
k_hyd_PU = 0.05125, J_hyd_max = 0.037625, K_Glc_PU = 8.7, nhyd_PU = 2.7;

G_max_VDDR = 3, G_max_bkc= 23;
tau_d_CaCl = .03, tau_d_NSCC = .35, G_max_NSCC = 12.15;

Na_o = 237, K_o = 7, 
K_i = 120, 
Na_K_Perm = 1.056075;

eta = .015;
Vm_IP3 = 3.33e-5, k4 = .0005, Pmv = 1.33e-5, kv = -58;

sigmaicc = .5, Amicc = 100, Cmicc = .01;
p_IP3SM = 5.3129e-04; P_gjSM = 1.4982e-11;

z_Ca = 2;       %Ca ion valence
P_gjICC = 0*(1.8458025e-12/3.27) % original is 3.27 cm^3/s
p_IP3ICC = 0*0.2*40  %1/s

%%SMC parameters

Ach = 0.00001, Gcouple = 1.3;
F_sm = 96.486, R_sm = 8.3144;
Q10K_sm = 1.365;
K_o_sm = 5.9, Cl_o = 134;
Cm_SM = .077, Vol_SM = 3500;
Na_i_sm = 10, K_i_sm = 164;
I_stim_period = 20, delta_VICC = 59;
zetaf = .5, k_CO = .01
Gmax_Ltyp = 65, Jmax_CASR = 0.31705;
Gmax_LVA = 0.18, Gmax_bkc = 45.7, Gmax_bk = 0.0144;
Gmax_kr = 35, Gmax_Na = 3, Gmax_ka = 9;
E_NSCC_sm = -28, Gmax_NSCC = 50;
E_K_sm = -88.8196, E_Na_sm = 69.9194;
gcoup = 0.5;



%%Equations

F_RT = F / (R * Temp);
RT_F = (R*Temp) / F;

F_RT_sm = F_sm / (R_sm * Temp);
RT_F_sm = (R_sm*Temp) / F_sm;

Vol_P_cyto = Vol_ICC * P_cyto_ICC;
Vol_P_PU = Vol_ICC * P_PU;
Vol_P_mito = Vol_ICC * P_mito_PU;
VOL_P_ER = Vol_ICC * P_ER_PU;

JhdSS = J_hyd_max/(1+((K_Glc_PU/Glc_PU)^(nhyd_PU)));


%% initial conditions
%% Initial Concentrations & transmembrane potential Vm
for i=1:number
    beta_IP3(i) = 2.67e-5 - .01e-5*(i-1);  %changes frequency
end       

gj_array = zeros(number,1);
for i=1:number
%     gj_array_ip3(i,1) = 8 + (-1+2*rand(1,1))*Mu*8
%     gj_array(i,1) = .53 + (-1*rand(1,1))*Mu*.7
    gj_array(i,1) = .70 + (-1+2*rand(1,1))*0.2*.7
end

Vm_in = -67*ones(1,number);
Ca_i_in = 0.00000993087*ones(1,number);
% K_i_in(1:number) = 8.4;
d_Na_in = 0*ones(1,number);
f_Na_in = 1*ones(1,number);
d_kv_in = 0*ones(1,number);
f_kv_in = 1*ones(1,number);
d_ERG_in = 0*ones(1,number);
Ca_PU_in = 0.0000902*ones(1,number);
d_Ltyp_in = 0*ones(1,number);
f_Ltyp_in = 1*ones(1,number);
fCa_Ltp_in = 1*ones(1,number);
ADP_i_in = 0.0077282*ones(1,number);
ADP_m_in = 2.60093454*ones(1,number);
h_PU_in = 0.9397*ones(1,number);
Ca_ER_PU_in = 0.007299*ones(1,number);
deltaPsi_in = 164.000044*ones(1,number);
d_VDDR_in = 0*ones(1,number);
f_VDDR_in = 1*ones(1,number);
d_CaCl_in = 0*ones(1,number);
d_NSCC_in = 0*ones(1,number);
NADH_m_in = 0.101476*ones(1,number);
Ca_m_PU_in = 0.000136*ones(1,number);
IP3_PU_in = 0.00064*ones(1,number);

Vm_sm_in = -69.75*ones(1,number);
Ca_i_sm_in = 0.00008*ones(1,number);
d_Ltyps_in = 0*ones(1,number);
f_Ltyps_in = 0.95*ones(1,number);
f_Ca_Ltyp_in = 1.0*ones(1,number);
d_LVA_in = .02*ones(1,number);
f_LVA_in = .99*ones(1,number);
xr1_in = 0*ones(1,number);
xr2_in = .82*ones(1,number);
m_Na_in = 0.005*ones(1,number);
h_Na_in = 0.05787*ones(1,number);
xa1_in = 0.00414*ones(1,number);
xa2_in = 0.72*ones(1,number);
m_NSCC_in = 0*ones(1,number);
 


%% Initial conditions vector
yinit= [Vm_in Ca_i_in d_Na_in f_Na_in d_kv_in f_kv_in d_ERG_in Ca_PU_in d_Ltyp_in f_Ltyp_in fCa_Ltp_in ADP_i_in ADP_m_in...
    h_PU_in Ca_ER_PU_in deltaPsi_in d_VDDR_in f_VDDR_in d_CaCl_in d_NSCC_in NADH_m_in Ca_m_PU_in IP3_PU_in...
    Vm_sm_in Ca_i_sm_in d_Ltyps_in f_Ltyps_in f_Ca_Ltyp_in d_LVA_in f_LVA_in xr1_in xr2_in m_Na_in h_Na_in xa1_in xa2_in m_NSCC_in];


%% Solve
%Default NE and NO stimulation: NEstimulation = 1e-3;   NOstimulation = 1e-3;   k_leak = 1;
%ameters for Oscillations:   NEstimulation = 0.4e-3; NOstimulation = 0;      k_leak = 5;
%k_pG = 0 to 0.1e-3 [1/ms] - receptor phosphorylation rate, keep k_pG = 0 to eliminate desensitization to NE;
% tic
[T,Y] = ode15s(@Eqns,tspan,yinit);
% toc 

%% Compute state variables
ii = 1;
Vm1 = Y(:,ii:number);
Ca_i1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 21;
% K_i1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_Na1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% f_Na1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_kv1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% f_kv1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_ERG1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% Ca_PU1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_Ltyp1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% f_Ltyp1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% fCa_Ltp1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% ADP_i1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% ADP_m1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% h_PU1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% Ca_ER_PU1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% deltaPsi1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_VDDR1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% f_VDDR1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_CaCl1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_NSCC1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% NADH_m1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% Ca_m_PU1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
IP3_PU1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;

Vm_sm1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
Ca_i_sm1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_Ltyps1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% f_Ltyps1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% f_Ca_Ltyp1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% d_LV1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% f_LVA1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% xr11 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% xr21 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% m_Na1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% h_Na1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% xa11 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% xa21 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% m_NSCC1 = Y(:,ii*number+1:(ii+1)*number); ii = ii + 1;

%% Compute non state variables (currents and reversal potentials)
% Is_and_Es = zeros(length(T),non_state_num);      % Currents and Reversal potentials
% 
% for jj = 1:length(T)
%     Eqns(T(jj),Y(jj,:));
%     Is_and_Es(jj,:) = non_state_vars;
% end
% 
% ii = 1;
% 
% I_Na    = Is_and_Es(:,ii:number);           %[mV] Reversal potentials
% % I_ktot  = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% % I_Catot = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% % I_CaCl    = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_gjICC = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% Itot    = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% save('267_42_7CaIP340X_1_both')



%%
    function dy = Eqns(t,y)

               
        %% Define state variables
        y = y(:);
        ii = 1;
        Vm	 	= y(ii:number)';
        Ca_i 	= y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_Na 	= y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        f_Na 	= y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_kv 	= y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        f_kv 	= y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_ERG = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        Ca_PU = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_Ltyp = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        f_Ltyp = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        fCa_Ltp = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        ADP_i = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        ADP_m = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        h_PU 	= y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        Ca_ER_PU 	= y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        deltaPsi 	= y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_VDDR = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        f_VDDR = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_CaCl = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_NSCC = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        NADH_m = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        Ca_m_PU = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        IP3_PU = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        
        Vm_sm = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        Ca_i_sm = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_Ltyps = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        f_Ltyps = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        f_Ca_Ltyp = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        d_LVA = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        f_LVA = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        xr1 = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        xr2 = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        m_Na = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        h_Na = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        xa1 = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        xa2 = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
        m_NSCC = y(ii*number+1:(ii+1)*number)'; ii = ii + 1;
      
        
        dy = zeros(size(y));
        state_num = length(y);
        
        
        
        %% Gap junctional currents and fluxes
        
        
        J_IP3 = zeros(number);
        J_elec = zeros(number);
        I_gjCa = zeros(number);
        
        
        for ICC_i = 1:number
            for ICC_j = 1:number
                if Index(ICC_i,ICC_j) ~= 0
                    % J_elec(ICC_i,ICC_j) = 0.7*(Vm(ICC_i)-Vm(ICC_j));
                    % J_elec(ICC_i,ICC_j) = (.7 + (-1+2*rand(1,1))*Mu*.7)*(Vm(ICC_i)-Vm(ICC_j));
                    J_elec(ICC_i,ICC_j) = gj_array(ICC_i,1)*(Vm(ICC_i)-Vm(ICC_j));
                    J_IP3(ICC_i,ICC_j) = p_IP3ICC*(IP3_PU(ICC_i)-IP3_PU(ICC_j));
                    if ((Vm(ICC_i) - Vm(ICC_j)) == 0)
                        I_gjCa(ICC_i,ICC_j)=1D12*P_gjICC*z_Ca*F*(Ca_i(ICC_i) - Ca_i(ICC_j));	%[pA]
                    else
                        I_gjCa(ICC_i,ICC_j)=1D12*P_gjICC *(Vm(ICC_i)-Vm(ICC_j))*(z_Ca^2)*F/RT_F*(Ca_i(ICC_j) - Ca_i(ICC_i)*exp((Vm(ICC_i)-Vm(ICC_j))*z_Ca/RT_F))/(1-exp((Vm(ICC_i)-Vm(ICC_j))*z_Ca/RT_F)); %[pA]
                    end
                end
            end
        end
        
        
        J_IP3ICC = sum(J_IP3,2)';
        I_gjICCIP3 =  1e12.*F.*J_IP3ICC.*Vol_ICC;
        
        I_gjICCelec = sum(J_elec,2)';
        I_gjICCCa = sum(I_gjCa,2)';
        
        I_gjICC =  I_gjICCelec + I_gjICCIP3 + I_gjICCCa;
        
        d_inf_Na  = 1./(1+exp((Vm +47)./(-4.8)));
        d_Naprime = (d_inf_Na -d_Na )./(.003.*Q10Na.^((Temp-T_exp)./10));
        f_inf_Na  = 1./(1+exp((Vm +78)./(7)));
        f_Naprime = (f_inf_Na -f_Na )./(.0016.*Q10Na.^((Temp-T_exp)./10));
        I_Na  = G_max_Na.*d_Na.*f_Na.*(Vm  - E_Na);
        
        d_inf_kv  = 1./(1+exp((Vm +25)./(-7.7)));
        d_kvprime = (d_inf_kv -d_kv )./(.005*Q10K^((Temp-T_exp)./10));
        f_inf_kv  = .5 + (.5./(1+exp((Vm +44.8)./(4.4))));
        f_kvprime = (f_inf_kv -f_kv )./(.005*Q10K^((Temp-T_exp)./10));
        I_kv  = G_max_kv.*d_kv.*f_kv.*(Vm  - E_K);

        dinf_ERG  = .2 +  (.8./(1+exp((Vm +20)./(-1.8))));
        d_ERGprime = (dinf_ERG -d_ERG )./(.003*Q10K^((Temp-T_exp)./10));
        I_ERG  = G_max_ERG.*d_ERG.*(Vm  - E_K);

        I_bk  = G_max_bk.*(Vm  - E_K);
        
        d_inf_L  = 1./(1+exp((Vm +17)./ (- 4.3)));
        d_Ltypprime = (d_inf_L -d_Ltyp )./(.001*Q10Ca^((Temp-T_exp)./10));
        f_inf_L  = 1./(1+exp((Vm +43)./8.9));
        f_Ltypprime =  (f_inf_L  - f_Ltyp )./(.086*Q10Ca^((Temp-T_exp)./10));
        fCa_infL  = 1 - (1./(1+exp((Ca_i  - 0.0001 - 0.000214)./ (- 1.310e-5))));
        fCa_Ltpprime = (fCa_infL -fCa_Ltp )./(.002*Q10Ca^((Temp-T_exp)./10));
        E_Ca_L  = 0.5*RT_F*log(Ca_o./Ca_i );
        I_Ltyp  = G_max_Ltyp.*d_Ltyp .*f_Ltyp .*fCa_Ltp .*(Vm  - E_Ca_L );


        
        J_leak  =  Jmax_leaK*(Ca_PU  - Ca_i );
        J_PMCA  = ( J_max_PMCA)./(1+(0.000298./Ca_i ));

        J_NaCa  = ( Jmax_NaCa*exp( b_PU*F_RT*(deltaPsi  - deltaPsi_s)))./( (1 +((K_Na_PU./Na_i).^n_PU)).*(1+(K_Ca_PU./Ca_m_PU )));

        PMF_PU  = deltaPsi  -  (2.303*RT_F*deltapH_PU);
        J_Hleak  = g_H_PU*PMF_PU;
        
       
        NAD_m_PU  = tot_NAD_m - NADH_m;
        A_res_PU  = RT_F.*log(( Kres_PU.*(NADH_m.^.5))./(NAD_m_PU.^.5));
        J_Hres  = ( rho_res_PU.*3.966.*( ra_PU.*10^ ( 6.*deltapH_PU).*exp( F_RT.*A_res_PU ))+ (rb_PU.*10^(6.*deltapH_PU))+  ((- 1).*(ra_PU+rb_PU).*exp(g_PU.*F_RT.*deltaPsi .*6))) ./ (((1+ r1_PU.*exp( F_RT.*A_res_PU )).*exp(6.*F_RT.*deltaPsi_B))+ ((r2_PU+ r3_PU.*exp(F_RT.*A_res_PU )).*exp(g_PU.*6.*F_RT.*deltaPsi )));

        f_PDHa_PU  = 1./(1 + (u2_PU.*(1+(u1_PU./((1+(Ca_m_PU ./KCa_PDH_PU)).^ 2)))));

        ATP_i  = tot_ANP_i - ADP_i;
        ATP_m  = tot_ANP_m - ADP_m; 

        J_glytot  = ( beta_max*(1+ beta1_PU*Glc_PU)*beta2_PU*Glc_PU*ATP_i )./(1+ (beta3_PU*ATP_i ) + ((1+ beta4_PU*ATP_i )*beta5_PU*Glc_PU) + ((1+ beta6_PU*ATP_i )*beta7_PU*Glc_PU));
        J_red_PU  = J_red_bas+ 6.39440.*f_PDHa_PU.*J_glytot; 
 
        
        J_hyd_PU  = k_hyd_PU*ATP_i +JhdSS;
        J_PGly  = 0.15*J_glytot;
        
        ADP_mfr  = 0.8*ADP_m; 
        ADP_ifr  = 0.3*ADP_i;

        ADP3_m  = 0.45*ADP_mfr; 
        ATP4_m  = 0.05*ATP_m;
        ADP3_i  = 0.45*ADP_ifr;  
        ATP4_i  =  0.05*ATP_i; 
        
        J_ANT_PU  = ( Jmax_ANT.*(1 -  ((( ATP4_i .*ADP3_m )./( ADP3_i .*ATP4_m )).*exp(  - F_RT.*deltaPsi ))))./( (1+ ((ATP4_i ./ADP3_i ).*exp(-frac_PU.*F_RT.*deltaPsi ))).*(1+(ADP3_m ./ATP4_m ))); 
        J_pTCA  = (J_red_bas./3) + (0.84.*f_PDHa_PU .*J_glytot );
        
        A_F1_PU  = RT_F*log(( K_F1_PU*ATP_m )./( ADP_mfr *Pi_m_PU));
        J_pF1_PU  = ( -rho_F1_PU.*( ((pa_PU.*10^  (3.*deltapH_PU))+ (pc1_PU.*exp( 3.*F_RT.*deltaPsi_B))).*exp( F_RT.*A_F1_PU ) - (pa_PU.*exp( 3.*F_RT.*deltaPsi ))+ (pc2_PU.*exp( F_RT.*A_F1_PU ).*exp( 3.*F_RT.*deltaPsi )))) ./ (((1+ p1_PU.*exp( F_RT.*A_F1_PU )).*exp( 3.*F_RT.*deltaPsi_B))+ ((p2_PU+ p3_PU.*exp( F_RT.*A_F1_PU )).*exp( 3.*F_RT.*deltaPsi )));
 
        
        ADP_iprime = ((-J_ANT_PU *Vol_P_mito)./Vol_P_cyto) + J_hyd_PU  - J_PGly; 
        ADP_mprime = J_ANT_PU  - J_pTCA  - J_pF1_PU; 

        J_o_PU  = ( rho_res_PU.*0.5.*(((( ra_PU.*10.^(6.*deltapH_PU))+ (rc1_PU.*exp( 6.*deltaPsi_B.*F_RT))).*exp( A_res_PU .*F_RT))  - (ra_PU.*exp( g_PU.*6.*F_RT.*deltaPsi ))+ (rc2_PU.*exp( F_RT.*A_res_PU ).*exp( F_RT.*deltaPsi .*6.*g_PU)))) ./ (((1+ r1_PU.*exp( F_RT.*A_res_PU )).*exp( F_RT.*deltaPsi_B.*6))+ ((r2_PU+ r3_PU.*exp( F_RT.*A_res_PU )).*exp( F_RT.*deltaPsi .*g_PU.*6)));
        NADH_mprime = J_red_PU  - J_o_PU;
 
        J_HF1_PU  = (  - rho_F1_PU.*3.*(( pa_PU.*10.^(3.*deltapH_PU).*exp( F_RT.*A_F1_PU ))+ (pb_PU.*10.^(3.*deltapH_PU)) - ((pa_PU+pb_PU).*exp( 3.*F_RT.*deltaPsi )))) ./ (((1+ p1_PU.*exp( F_RT.*A_F1_PU )).*exp(3.*F_RT.*deltaPsi_B))+ ((p2_PU+ p3_PU.*exp( F_RT.*A_F1_PU )).*exp( 3.*F_RT.*deltaPsi )));
        MWC_PU  = ( (( conc_PU.*Ca_PU )./K_trans_PU).*((1+(Ca_PU ./K_trans_PU)).^ 3))./(((1+(Ca_PU ./K_trans_PU)).^ 4)+(L_PU./((1+(Ca_PU ./K_act_PU)).^ Na_PU)));
        J_uni_PU  = ( Jmax_uni.*(MWC_PU  -  (Ca_m_PU .*exp(  - 2.*F_RT.*(deltaPsi  - deltaPsi_s)))).*2.*F_RT.*(deltaPsi  - deltaPsi_s))./(1 - exp(  - 2.*F_RT.*(deltaPsi  - deltaPsi_s)));


        J_ERout  =  ( Jmax_IP3.*((IP3_PU ./(IP3_PU +d_IP3_PU)).^ 3).*((Ca_PU ./(Ca_PU +d_ACT_PU)).^ 3).*(h_PU .^ 3)+J_ERleak).*(Ca_ER_PU  - Ca_PU );
        J_serca  = ( Jmax_serca.*(Ca_PU .^ 2))./((K_serca_PU.^ 2)+(Ca_PU .^ 2));

        h_PUprime = (d_INH_PU -  (h_PU .*(Ca_PU +d_INH_PU)))./tauh_PU;
        Ca_ER_PUprime = fe_PU.*(J_serca  - J_ERout ); 
        Ca_m_PUprime = fm_PU.*(J_uni_PU  - J_NaCa );


        d_bkc  = 1./(1+exp((Vm ./ (- 17)) -  (2.*log(Ca_i ./0.001))));
        I_bkc  = (G_max_bkc+1.1.*(Temp-T_exp)).*d_bkc .*(Vm  - E_K);

        DInfCaCl  = 1./(1+(.00014./Ca_i ).^3);
        d_CaClprime = ((DInfCaCl  - d_CaCl )./tau_d_CaCl);
        I_CaCl  = G_max_CaCl.*d_CaCl .*(Vm  - E_Cl);

        DInfVDDR  = 1./(1+exp((Vm +26)./ (-6)));
        d_VDDRprime =  (DInfVDDR  - d_VDDR )./(.006.*Q10Ca.^((Temp-T_exp)./10));
        FInfVDDR  = 1./(1+exp((Vm +66)./6));
        f_VDDRprime = (FInfVDDR  - f_VDDR )./(.04.*Q10Ca.^((Temp-T_exp)./10));
        E_CaVDDR = 0.5.*RT_F.*log(Ca_o./Ca_i );
        I_VDDR  = G_max_VDDR.*d_VDDR .*f_VDDR .*(Vm  - E_CaVDDR );

        DInfNSCC  = 1./(1+((7.45000e-05./Ca_PU ).^(- 85)));
        d_NSCCprime = (DInfNSCC  - d_NSCC )./tau_d_NSCC;
        E_NSCC = RT_F.*log((K_o+ (Na_o.*Na_K_Perm))./(K_i+ (Na_i.*Na_K_Perm)));
        I_NSCC  =  G_max_NSCC.*d_NSCC .*(Vm  - E_NSCC);

       
        
       
        deltaPsiprime = (((  -F.*Vol_P_mito.*1e6)./Cmito_PU).*(J_Hleak -J_Hres +J_ANT_PU +J_HF1_PU + (2.*J_uni_PU )));

        Ca_iprime = fc_ICC.*(((-I_Ltyp-I_VDDR-I_gjICCCa)./( 2e12.*F.*Vol_ICC .* P_cyto_ICC)) + J_leak - J_PMCA );
%         K_iprime = fc_ICC.*(((-I_kv-I_ERG-I_bk-I_bkc-I_gjICCK)./( 1e12.*F.*Vol_ICC .* P_cyto_ICC)));
        Ca_PUprime = fc_ICC.*(((( J_NaCa  - J_uni_PU ).*Vol_P_mito)./Vol_P_PU) +(((J_ERout  - J_serca ).*VOL_P_ER)./Vol_P_PU) + ((  -J_leak .*Vol_P_cyto)./Vol_P_PU));

        IP3_PUprime = beta_IP3 - eta.*IP3_PU  - ((Vm_IP3.*(IP3_PU.^4))./(k4^4+IP3_PU.^4)) + Pmv.*(1 - (Vm.^8./(kv^8 + Vm.^8)))- J_IP3ICC;
        
                   
                 
        I_ktot = I_kv + I_ERG + I_bk + I_bkc;
        I_Catot = I_Ltyp+I_VDDR;
        Itot = I_Na+I_Ltyp+I_VDDR+I_ERG+I_bk+I_kv+I_bkc+I_CaCl+I_NSCC+J_PMCA.*2e12.*F.*Vol_P_cyto;
        
        
        
 
 %% Smooth muscle equations
 
dinf_Ltyps = 1./(1+exp((Vm_sm+17)./(-4.3)));
tau_d_Lts = .00047.* ((Q10Ca).^((Temp-T_exp)./10));
d_Ltypsprime = 1.*((dinf_Ltyps - d_Ltyps)./ tau_d_Lts);
 
finf_Ltyps = 1./(1+exp((Vm_sm+43)./(8.9)));
tau_f_Lts = .086 .* ((Q10Ca).^((Temp-T_exp)./10));
f_Ltypsprime = 1.*((finf_Ltyps - f_Ltyps) ./ tau_f_Lts);
 
finf_Ca_Lt = 1./(1+exp((Ca_i_sm-0.00008999-0.000214)./(-0.0000131)));
tau_fCa_Lt = .002.* ((Q10Ca).^((Temp-T_exp)./10));
f_Ca_Ltypprime = 1.*((finf_Ca_Lt - f_Ca_Ltyp) ./ tau_fCa_Lt);
 
dinf_LVA = 1./(1+exp((Vm_sm+27.5)./(-10.9)));
tau_d_LVA = .003 .* ((Q10Ca).^((Temp-T_exp)./10));
d_LVAprime = 1.*((dinf_LVA - d_LVA) ./ tau_d_LVA);
 
finf_LVA = 1./(1+exp((Vm_sm+15.8)./(7)));
tau_f_LVA = .00758 .* ((Q10Ca)^((Temp-T_exp)./10)) .* exp(Vm_sm.*0.00817);
f_LVAprime = 1.*((finf_LVA - f_LVA) ./ tau_f_LVA); 
 
xr1_inf = 1./(1+exp((Vm_sm+27)./(-5.0)));
tau_xr1 = .080 .* ((Q10K_sm).^((Temp-T_exp)./10));
xr1prime = 1.*((xr1_inf - xr1) ./ tau_xr1);
 
xr2_inf = .8./(1+exp((Vm_sm+58)./(10.0)));
tau_xr2 = ((Q10K_sm).^((Temp-T_exp)./10)) .* (-.707 + 1.481.*exp((Vm_sm+36)./92));
xr2prime = 1.*((xr2_inf - xr2) ./ tau_xr2);
 
xa1_inf = 1./(1+exp((Vm_sm+26.5)./(-7.9)));
tau_xa1 = ((Q10K_sm).^((Temp-T_exp)./10)) .* (.0318 + .175.*exp((-.5).*(((Vm_sm+44.4)./22.3).^2)));
xa1prime = 1.*((xa1_inf - xa1) ./ tau_xa1);
 
xa2_inf = .1 +  (.9./(1+exp((Vm_sm+65)./(6.2))));
tau_xa2 = .090 .* ((Q10K_sm).^((Temp-T_exp)./10));
xa2prime = 1.*((xa2_inf - xa2) ./ tau_xa2);
 
 
minf_Na = 1./(1+exp((Vm_sm+47)./(-4.8)));
tau_m_Na = ((Q10Na).^((Temp-T_exp)./10)) .* ((-.000017.*Vm_sm) + .00044);
m_Naprime = 1.*((minf_Na- m_Na) ./ tau_m_Na);
 
hinf_Na = 1./(1+exp((Vm_sm+78)./(3.0)));
tau_h_Na = ((Q10Na).^((Temp-T_exp)./10)) * ((-.00025.*Vm_sm) + .0055);
h_Naprime = 1.*((hinf_Na - h_Na) ./ tau_h_Na);
 
minf_NSCC = 1./(1+exp((Vm_sm+25)./(-20)));
tau_m_NSCC = .150./(1+exp((Vm_sm+66)./(-26)));
m_NSCCprime = 1.*((minf_NSCC- m_NSCC) ./ tau_m_NSCC);
 
E_Ca = .5 .* RT_F_sm .* log(Ca_o./Ca_i_sm);
I_Ltyp_sm = Gmax_Ltyp .* d_Ltyps .* f_Ltyps .* f_Ca_Ltyp .* (Vm_sm - E_Ca);
I_LVA = Gmax_LVA .* f_LVA .* d_LVA .* (Vm_sm - E_Ca);
 
I_Na_sm = Gmax_Na .* m_Na .* h_Na .* (Vm_sm-E_Na_sm);

% CO_conc = .1 + .4./(1+exp(-(zeta - zetaf)./k_CO));
% fCO = 2.475.*CO_conc - .2375;
I_kr = Gmax_kr .* xr1 .* xr2 .* (Vm_sm-E_K_sm);

I_ka = Gmax_ka .* xa1 .* xa2 .* (Vm_sm-E_K_sm);
 
d_bkc_sm = 1./(1 + exp((-Vm_sm./17)-2.*log(Ca_i_sm./.001)));
I_bkc_sm = (Gmax_bkc + 1.1.*(Temp-T_exp)) .* d_bkc_sm .* (Vm_sm-E_K_sm);
 
fCa_NSCC = 1./(1+((Ca_i_sm./.0002).^(-4)));
rach_NSCC = 1./(1 + (.01./Ach));
I_NSCC_sm = Gmax_NSCC .* m_NSCC .* fCa_NSCC .* rach_NSCC .* (Vm_sm - E_NSCC_sm);
 
I_bk_sm = Gmax_bk .* (Vm_sm - E_K_sm);
 
J_CASR = 1000 .*Jmax_CASR .* ((Ca_i_sm).^1.34);
Ca_i_smprime = 1.*((-I_Ltyp_sm-I_LVA)./(2.*F_sm.*.001.*Vol_SM)) - J_CASR;
 
        
Itot_sm =  I_LVA+I_Ltyp_sm+I_Na_sm+I_kr+I_ka+I_bk_sm+I_bkc_sm+I_NSCC_sm;
I_coup= gcoup.*(Vm-Vm_sm); %%Mu is replacement of gcoup

Vmprime = (-1/Cmicc).*(Itot + I_gjICC );
Vm_smprime = (-1./Cm_SM).*(Itot_sm - I_coup);
  
       
        
       
             
        %% Differential Equations
        ii = 1;
        dy(ii:number)  = Vmprime(1:number);
        dy(ii*number+1:(ii+1)*number)  = Ca_iprime(1:number); ii = ii + 1;
%         dy(ii*number+1:(ii+1)*number)  = K_iprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number)  = d_Naprime(1:number) ; ii = ii + 1;
        dy(ii*number+1:(ii+1)*number)  = f_Naprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number)  = d_kvprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number)  = f_kvprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number)  = d_ERGprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = Ca_PUprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number)  = d_Ltypprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = f_Ltypprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = fCa_Ltpprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = ADP_iprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = ADP_mprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = h_PUprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = Ca_ER_PUprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = deltaPsiprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = d_VDDRprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = f_VDDRprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = d_CaClprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = d_NSCCprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = NADH_mprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = Ca_m_PUprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = IP3_PUprime(1:number); ii = ii + 1;
        
        dy(ii*number+1:(ii+1)*number) = Vm_smprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = Ca_i_smprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = d_Ltypsprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = f_Ltypsprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = f_Ca_Ltypprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = d_LVAprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = f_LVAprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = xr1prime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = xr2prime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = m_Naprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = h_Naprime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = xa1prime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = xa2prime(1:number); ii = ii + 1;
        dy(ii*number+1:(ii+1)*number) = m_NSCCprime(1:number); ii = ii + 1;
       
        
            
       
        %% Non state variables
%         non_state_vars = [I_Na I_gjICC Itot];
%         non_state_num = length(non_state_vars);
        
    end
toc
end



