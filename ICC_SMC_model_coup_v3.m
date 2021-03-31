function ICC_model
tic
clc, close all
warning('off')

%% Define simulation time
TMAX  = 600;   %[ms] simulation time
tspan = [0 TMAX];

%% Total number of ICCs
number = 11;

%% Define connectivity (index) matrix
parent = 0:number-1;
n = length(parent);
m = find(parent);
Index = sparse (parent(m), m, 1, n, n);
Index = full(Index);
Index = Index+triu(Index,1).';

%% make graph object
Adj = Index;
Gobj = graph(Adj);

%% find neighbors
nnod = numnodes(Gobj);

edgelist = Gobj.Edges{:,1};
neighbors = zeros(nnod, max(Gobj.degree));
for i = 1:nnod
    neighborNodes = Gobj.neighbors(i);
    neighborNodes = neighborNodes';
    neighbors(i,1:numel(neighborNodes)) = neighborNodes;
end

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
% p_IP3SM = 5.3129e-04; P_gjSM = 1.4982e-11;

z_Ca = 2;       %Ca ion valence
P_gjICC = 1.8458025e-12/3.27 % original is 3.27 cm^3/s
p_IP3ICC = 8.0  %1/s
p_CaICC = 0.1

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
%     zeta(i) = 0 + .2*(i-1);
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
tic
disp('------------ Finding jacobian matrix pattern ------------')

Adj = Index;
A = Adj + eye(size(Adj));
rep = numel(yinit)/number;
A = repmat(A,[rep,rep]);
JPattern = A;
JPattern = sparse(JPattern);
toc

opts=odeset('JPattern',JPattern);
tic
[T,Y] = ode15s(@Eqns,tspan,yinit,opts);
toc

%% Compute state variables
ii = 1;
Vm_icc1 = Y(:,ii:number);
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

% idx = find(T>=TMAX-40)
% 
% Is_and_Es = zeros(length(idx),non_state_num);      % Currents and Reversal potentials
% for jj = idx(1):idx(end)
%     Eqns(jj,Y(jj,:));
%     Is_and_Es(jj,:) = non_state_vars;
% end
% Is_and_Es = Is_and_Es(length(Is_and_Es)-length(idx):end,:);

% Is_and_Es = zeros(length(T),non_state_num); 
% for jj = 1:length(T)
%     Eqns(T(jj),Y(jj,:));
%     Is_and_Es(jj,:) = non_state_vars;
% end
% 
% ii = 1;
% 
% I_Na    = Is_and_Es(:,ii:number);           %[mV] Reversal potentials
% I_Ltyp  = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_VDDR = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_CaCl    = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_kv    = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_ERG    = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_bk    = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_bkc    = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_NSCC    = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_gjICC = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% 
% I_Ltyp_sm = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_LVA = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_Na_sm = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_kr = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_ka = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_bkc_sm = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_bk_sm = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_NSCC_sm = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% I_coup = Is_and_Es(:,ii*number+1:(ii+1)*number); ii = ii + 1;
% 
% 
% I_NSCC_in = zeros(length(I_NSCC),number);
% I_NSCC_out = zeros(length(I_NSCC),number);
% for i = 1:length(I_NSCC)
%     for j=1:number
%         if I_NSCC(i,j)>=0
%             I_NSCC_in(i,j) = 0;
%             I_NSCC_out(i,j) = I_NSCC(i,j);
%         else
%             I_NSCC_in(i,j) = I_NSCC(i,j);
%             I_NSCC_out(i,j) = 0;
%         end
%     end
% end
clearvars Y
save('abc')

%% Plot structure, Vm and ionic concentrations vs. time
% tplot = T;
% T = T_new{1,1};
N = 1:number;

figure, plot(T,Vm_icc1(:,number),'linewidth',1)
hold on, plot(T,Vm_sm1(:,number),'linewidth',1),
ylim([-70 -20]);
legend('ICC cell','SM cell')
xlabel('Time (sec)'); ylabel('Membrane potential (mV)')
% set(gcf,'Position',[100 100 1500 500])
% 
% figure
% hold on, plot(T,Vm_sm1(:,number),'linewidth',1),
% ylim([-70 -40])
% xlabel('Time (sec)'); ylabel('Membrane potential (mV)')
% set(gcf,'Position',[100 100 1500 500])

figure
subplot(2,1,1)
plot(T,Vm_icc1(:,1),T,Vm_icc1(:,number))
legend('1st ICC cell','Last ICC cell')
subplot(2,1,2)
plot(T,Vm_sm1(:,1),T,Vm_sm1(:,number))
legend('1st SM cell','Last SM cell')
% plot(T,Ca_i1(:,1),'r',T,Ca_i_sm1(:,1),'b')
% plotmatrix(Vm_icc1(:,1),Vm_icc1(:,number),'ob')

%Blowup
a = find(T> max(T)-120); %%Find the index of last 120 sec of data
b = T(T> max(T)-120); %% Find exact time points in the selected region for plotting data
figure(3)
plot(b,Vm_sm1(a,1),b,Vm_sm1(a,number),'linewidth',1)
box off
% set(gca,'xticklabel',{[]});set(gca,'yticklabel',{[]})
% set(gca,'XMinorTick','off','YMinorTick','off')
legend('1st SM cell','Last SM cell')
xlabel('Time (sec)');ylabel('Membrane potential (mV)')



% % Currentscape diagram
% indx = 1:length(T);
%         
% I_in = abs(I_Ltyp(indx,1)) + abs(I_VDDR(indx,1)) + abs(I_CaCl(indx,1)) + abs(I_Na(indx,1)) + abs(I_NSCC_in(indx,1));
% I_out = -(I_kv(indx,1)) -(I_ERG(indx,1)) -(I_bk(indx,1)) -(I_bkc(indx,1)) -(I_NSCC_out(indx,1));
% I_tot_in = [abs(I_Ltyp(indx,1))./I_in(indx,1) abs(I_VDDR(indx,1))./I_in(indx,1)...
%         abs(I_CaCl(indx,1))./I_in(indx,1) abs(I_Na(indx,1))./I_in(indx,1) abs(I_NSCC_in(indx,1))./I_in(indx,1) ];
% I_tot_out = [I_kv(indx,1)./I_out(indx,1) I_ERG(indx,1)./I_out(indx,1)...
%         I_bk(indx,1)./I_out(indx,1) I_bkc(indx,1)./I_out(indx,1) I_NSCC_out(indx,1)./I_out(indx,1) ];
% 
% % https://www.mathworks.com/matlabcentral/fileexchange/32884-plot-groups-of-stacked-bars
% stackData_in = I_tot_in;
% stackData_out = I_tot_out;
% T = T';
% NumStacksPerGroup1 = size(stackData_in, 2);
% NumStacksPerGroup2 = size(stackData_out, 2);
% 
% 
% h=figure
%     hold on; 
% for i=1:length(T)-1
%     groupOffset = T(i+1) - T(i);
%     % Offset the group draw positions:
%     groupDrawPos =  T(i) + 0.5*(T(i+1)-T(i));
%     
%     h_in(i,:) = bar([stackData_in(i,:);zeros(1,NumStacksPerGroup1)], 'stacked','Edgecolor','none','FaceColor','flat');
%     h_out(i,:) = bar([stackData_out(i,:);zeros(1,NumStacksPerGroup2)], 'stacked','Edgecolor','none','FaceColor','flat');
%     h_in(i,1).CData = [1 0 0];
%     h_in(i,2).CData = [0 1 0];
%     h_in(i,3).CData = [0 0 1];
%     h_in(i,4).CData = [0 0 0];
%     h_in(i,5).CData = [.7 .7 .7];
%     h_out(i,1).CData = [1 1 0];
%     h_out(i,2).CData = [1 0 1];
%     h_out(i,3).CData = [0 1 1];
%     h_out(i,4).CData = [.3 .3 .3];
%     h_out(i,5).CData = [.7 .7 .7];
%     set(h_in(i,:),'BarWidth',groupOffset);
%     set(h_in(i,:),'XData',groupDrawPos);
%     set(h_out(i,:),'BarWidth',groupOffset);
%     set(h_out(i,:),'XData',groupDrawPos);
% end
% legend('L-type','VDDR','Cl','Na','NSCC','kv','ERG','bk','bkc'); box off
% % xlabel('Time (sec)');ylabel('% Contribution')
% xlim([0 20]); ylim([-1.2 1.2]); hold off;
% % set(gcf, 'Renderer', 'painters');
% savefig(h,'abc_compact.fig','compact')
% 
% I_in_sm = abs(I_Ltyp_sm(indx,1)) + abs(I_LVA(indx,1)) + abs(I_Na_sm(indx,1)) + abs(I_NSCC_sm(indx,1)) + abs(I_coup(indx,1));
% I_out_sm = -(I_ka(indx,1)) -(I_kr(indx,1)) -(I_bk_sm(indx,1)) -(I_bkc_sm(indx,1));
% I_tot_in_sm = [abs(I_Ltyp_sm(indx,1))./I_in_sm(indx,1) abs(I_LVA(indx,1))./I_in_sm(indx,1)...
%         abs(I_Na_sm(indx,1))./I_in_sm(indx,1) abs(I_NSCC_sm(indx,1))./I_in_sm(indx,1) abs(I_coup(indx,1))./I_in_sm(indx,1) ];
% I_tot_out_sm = [I_ka(indx,1)./I_out_sm(indx,1) I_kr(indx,1)./I_out_sm(indx,1)...
%         I_bk_sm(indx,1)./I_out_sm(indx,1) I_bkc_sm(indx,1)./I_out_sm(indx,1)];
% 
% % https://www.mathworks.com/matlabcentral/fileexchange/32884-plot-groups-of-stacked-bars
% stackData_in_sm = I_tot_in_sm;
% stackData_out_sm = I_tot_out_sm;
% 
% NumStacksPerGroup1_sm = size(stackData_in_sm, 2);
% NumStacksPerGroup2_sm = size(stackData_out_sm, 2);
% 
% 
% figure
%     hold on; 
% for i=1:length(T)-1
%     groupOffset_sm = T(i+1) - T(i);
%     % Offset the group draw positions:
%     groupDrawPos_sm =  T(i) + 0.5*(T(i+1)-T(i));
%     
%     h_in_sm(i,:) = bar([stackData_in_sm(i,:);zeros(1,NumStacksPerGroup1_sm)], 'stacked','Edgecolor','none','FaceColor','flat');
%     h_out_sm(i,:) = bar([stackData_out_sm(i,:);zeros(1,NumStacksPerGroup2_sm)], 'stacked','Edgecolor','none','FaceColor','flat');
%     h_in_sm(i,1).CData = [1 0 0];
%     h_in_sm(i,2).CData = [0 1 0];
%     h_in_sm(i,3).CData = [0 0 1];
%     h_in_sm(i,4).CData = [0 0 0];
%     h_in_sm(i,5).CData = [.7 .7 .7];
%     h_out_sm(i,1).CData = [1 1 0];
%     h_out_sm(i,2).CData = [1 0 1];
%     h_out_sm(i,3).CData = [0 1 1];
%     h_out_sm(i,4).CData = [.3 .3 .3];
% %     h_out(i,5).CData = [.7 .7 .7];
%     set(h_in_sm(i,:),'BarWidth',groupOffset_sm);
%     set(h_in_sm(i,:),'XData',groupDrawPos_sm);
%     set(h_out_sm(i,:),'BarWidth',groupOffset_sm);
%     set(h_out_sm(i,:),'XData',groupDrawPos_sm);
% end
% legend('L-type','LVA','Na','NSCC','Coup','ka','kr','bk','bkc'); box off
% % xlabel('Time (sec)');ylabel('% Contribution')
% xlim([0 20]); ylim([-1.2 1.2]); hold off;


figure 
subplot(7,1,1),plot(b,Vm_sm1(a,1),'linewidth',2);box off
for i=2:7
subplot(7,1,i),plot(b,Vm_sm1(a,(i-1)*7),'linewidth',2); box off
end

[pks1,locs1] = findpeaks(Vm_sm1(:,1),T(:,1),'MinPeakDistance',12);
[pks2,locs2] = findpeaks(Vm_sm1(:,number),T(:,1),'MinPeakDistance',12);

% https://www.mathworks.com/matlabcentral/answers/44227-finding-local-minimums-maximums-for-a-set-of-data
DataInv1 = max(Vm_sm1(:,1)) - Vm_sm1(:,1);
[valley1,MinIdx1] = findpeaks(DataInv1,T(:,1),'MinPeakDistance',12);
valley1 = max(Vm_sm1(:,1)) - valley1;

DataInv2 = max(Vm_sm1(:,number)) - Vm_sm1(:,number);
[valley2,MinIdx2] = findpeaks(DataInv2,T(:,1),'MinPeakDistance',12);
valley2 = max(Vm_sm1(:,number)) - valley2;

p_to_p1 = pks1(end-6:end) - valley1(end-6:end)
p_to_p2 = pks2(end-6:end) - valley2(end-6:end)


% locs2 = locs2(pks2>-51);
for i=1:length(locs2)
    phase_lag(i) = locs2(i) - locs1(i);
end
figure
plot(phase_lag(1:length(locs2)),'.','Markersize',30,...
    'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6]); box off
% plot(phase_lag(length(locs2)-6:length(locs2)),'linewidth',2);
xlabel('No of Cycle');ylabel('Lag (sec)')
ylim([0 220])
xlim([0 50])


figure
plot(length(pks2)-6:length(pks2),phase_lag(length(locs2)-6:length(locs2)),'.','Markersize',30,...
    'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6]);
% plot(1:length(pks2),phase_lag(1:length(locs2)),'.','Markersize',30,...
%     'MarkerEdgeColor','red', 'MarkerFaceColor',[1 .6 .6]);
ax = gca; box off;
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
xlabel('No of Cycle');ylabel('SWperiod (sec)')
% ylim([20.85 21.55])
% yticks([20.9 21.1 21.3 21.5])


for i=1:length(phase_lag)-1
    diff_phase_lag(i) = phase_lag(i+1) - phase_lag(i);
    period_first(i) = locs1(i+1) - locs1(i);
    period_last(i) = locs2(i+1) - locs2(i);
end
figure
plot(diff_phase_lag,'linewidth',1); ylim([-1 ceil(max(diff_phase_lag))]); box off
line([0,length(phase_lag)],[0,0],'Color','k') %%https://www.mathworks.com/matlabcentral/answers/306204-add-various-horizontal-lines-to-a-plot
xlabel('Cycle #');ylabel('Delta Lag (sec)')

figure
max_period = max(period_first(length(period_first)-6:length(period_first)));
plot(period_first(length(period_first)-6:length(period_first))/max_period,'*','Markersize',10);
hold on; plot(period_last(length(period_last)-6:length(period_last))/max_period,'+','Markersize',10)
xlabel('Cycle #');ylabel('Normalized Period (sec)'); ylim([0.5 1.5]); box off
legend('First SM cell','Last SM cell')

period_first_avg = mean(period_first(length(period_first)-6:length(period_first)));
period_last_avg = mean(period_last(length(period_first)-6:length(period_first)));
std_first_avg = std(period_first(length(period_first)-6:length(period_first)));
std_last_avg = std(period_last(length(period_first)-6:length(period_first)));
mean_freq_first = mean((1./period_first(length(period_first)-6:length(period_first)))*60);
mean_freq_last = mean((1./period_last(length(period_last)-6:length(period_last)))*60);
std_freq_first = std((1./period_first(length(period_first)-6:length(period_first)))*60);
std_freq_last = std((1./period_last(length(period_last)-6:length(period_last)))*60);
fprintf('Mean_frequency_1st_cell:%f\t',mean_freq_first);fprintf('Std_frequency_1st_cell:%f\n',std_freq_first);
fprintf('Mean_frequency_last_cell:%f\t',mean_freq_last);fprintf('Std_frequency_after:%f\n',std_freq_last);

%% Plot heat maps
figure
% number = 34
x = b;
y = 1:number
D(:,y) = Vm_sm1(a,y);
imagesc(x,y,D')
xlabel('Time(sec)');ylabel('Cell index')
% xticklabels({'1^s^t','7^t^h','14^t^h','21^s^t','28^t^h','35^t^h','42^n^d'})
colormap gray; colorbar
caxis([-64 -48])
% set(gcf,'Position',[100 100 700 500])

for i = 1:number
[pks,locs] = findpeaks(Vm_sm1(:,i),T(:,1),'MinPeakDistance',12);
locas{i,:} = locs;
end

len = size(locas{number,1},1)-1
for i=1:number
    locas{i,1} = locas{i,1}(1:len);
end

locs_new = zeros(number, len);
lags_new = zeros(number-1, len);

for i=1:number
    locs_new(i,:) = locas{i,1};
end

for i=1:number-1
    lags_new(i,:) = locs_new(i+1,:) - locs_new(i,:);
    lags_new_from_start(i,:) = locs_new(i+1,:) - locs_new(1,:);
end

% Lag diagram
figure
a1 = (1:len)';
b1=1:number-1;
% imagesc(a1,b1,lags_new)
imagesc(a1(end-20:end),b1,lags_new(:,end-20:end))
xlabel('Cycle #');ylabel('Cell Index')
colormap hot; colorbar
caxis([0 20])

figure
a1 = (1:len)';
b2=2:number;
imagesc(a1(end-20:end),b2,lags_new_from_start(:,end-20:end))
xlabel('Cycle #');ylabel('Cell Index')
colormap hot; colorbar
caxis([0 150])

% figure
% x=1:size(lags_new_from_start,1)
% plot(lags_new_from_start(1,end-20:end),'b')
% text(20,1,num2str(1))
% hold on
% % for i=1:size(lags_new_from_start,1)
% for i=2:6
%     plot(lags_new_from_start(7*i-7,end-20:end),'b')
%     text(20,max(lags_new_from_start(7*i-7,:)),num2str(7*i-7))
%     hold on
% end
% hold on; plot(lags_new_from_start(41,end-20:end),'b')
% text(20,1,num2str(41))

figure
for i=1:size(lags_new_from_start,1)
    plot(lags_new_from_start(i,end-20:end),'b')
    text(22,max(lags_new_from_start(i,end-20:end)),num2str(i))
    hold on
end
% ylim([0 220])

% Delta lag diagram
for i=1:len-1
    delta_lags_new(:,i) = lags_new(:,i+1) - lags_new(:,i);
end
figure
a2 = (1:len-1)';
% delta_lags_scaled = delta_lags_new(a1,b2)';
% imagesc(a2,b1,delta_lags_new)
imagesc(a2(end-20:end),b1,delta_lags_new(:,end-20:end))
xlabel('Cycle #');ylabel('Cell Index')
% title('Delta Lag')
colormap hot; colorbar
caxis([-1 4])

% for i=1:len-1
%     delta_lags_new_from_start(:,i) = lags_new_from_start(:,i+1) - lags_new_from_start(:,i);
% end
% figure
% a2 = (1:len-1)';
% % delta_lags_scaled = delta_lags_new(a1,b2)';
% % imagesc(a2,b2,delta_lags_new_from_start)
% imagesc(a2(end-20:end),b2,delta_lags_new_from_start(:,end-20:end))
% xlabel('Cycle #');ylabel('Cell Index')
% % title('Delta Lag')
% colormap hot; colorbar
% set(gcf,'Position',[100 100 700 500])

% Periodogram
for i=1:len-1
    period_new(:,i) = locs_new(:,i+1) - locs_new(:,i);
end

b3 = 1:number
figure
% imagesc(a1,b3,period_new)
imagesc(a1(end-20:end),b3,period_new(:,end-20:end))
xlim([len-20 len])
xlabel('Cycle #');ylabel('Cell Index')
% title('Period')
colormap hot; colorbar
caxis([17 21])
% set(gcf,'Position',[100 100 700 500])


%% Sensitivity analysis
% load('diff_phase_lag_267_28_3point5CaIP3_1_both.mat')
% a = delta_lag
% load('diff_phase_lag_267_28_7CaIP3_1_both.mat')
% b = delta_lag
% load('diff_phase_lag_267_28_14CaIP3_1_both.mat')
% c = delta_lag
% plot(a(length(a)-7:end))
% ylim([-.5 .5])
% hold on, plot(b(length(b)-7:end))
% hold on, plot(c(length(c)-7:end));box off
% xlabel('Cycle #');ylabel('Delta Lag (sec)');
% legend('0.5 nS', '1.0 nS', '2.0 nS')


% figure(5), plot(Vm_icc1(:,4)-Vm_icc1(:,5),I_gjICC(:,5))
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
        
        
        J_IP3 = zeros(number,1);
        J_elec = zeros(number,1);
        J_gjCa = zeros(number,1);
        
        for i = 1:number
            for kk = 1:nnz(neighbors(i,:))
                j = neighbors(i,kk);
                J_IP3(i) = J_IP3(i) + p_IP3ICC*(IP3_PU(j)-IP3_PU(i));
                J_elec(i) = J_elec(i) + 0.7*(Vm(i)-Vm(j));
                J_gjCa(i) = J_gjCa(i) + p_CaICC*(Ca_i(j)-Ca_i(i));
            end
        end
        
                    
        
        J_IP3ICC = sum(J_IP3,2)';
        I_gjICCIP3 =  1e12.*F.*J_IP3ICC.*Vol_ICC;
        
        I_gjICCelec = J_elec';
        J_gjICCCa = sum(J_gjCa,2)';
        
        I_gjICC =  I_gjICCelec;
        
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
        J_PMCA  = ( J_max_PMCA)./(1+(0.000298./Ca_i));

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

        Ca_iprime = fc_ICC.*(((-I_Ltyp-I_VDDR)./( 2e12.*F.*Vol_ICC .* P_cyto_ICC)) + J_leak - J_PMCA - J_gjICCCa);
%         K_iprime = fc_ICC.*(((-I_kv-I_ERG-I_bk-I_bkc-I_gjICCK)./( 1e12.*F.*Vol_ICC .* P_cyto_ICC)));
        Ca_PUprime = fc_ICC.*(((( J_NaCa  - J_uni_PU ).*Vol_P_mito)./Vol_P_PU) +(((J_ERout  - J_serca ).*VOL_P_ER)./Vol_P_PU) + ((  -J_leak .*Vol_P_cyto)./Vol_P_PU));

        IP3_PUprime = beta_IP3 - eta.*IP3_PU  - ((Vm_IP3.*(IP3_PU.^4))./(k4^4+IP3_PU.^4)) + Pmv.*(1 - (Vm.^8./(kv^8 + Vm.^8))) + J_IP3ICC;
        
                   
                 
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
I_coup= gcoup.*(Vm-Vm_sm);

Vmprime = (-1/Cmicc).*(Itot + I_gjICC);
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
%         non_state_vars = [I_Na I_Ltyp I_VDDR I_CaCl I_kv I_ERG I_bk I_bkc I_NSCC I_gjICC...
%                I_Ltyp_sm I_LVA I_Na_sm I_kr I_ka I_bkc_sm I_bk_sm I_NSCC_sm I_coup];
%         non_state_num = length(non_state_vars);
%         
    end
toc
end



