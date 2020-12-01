clear all; clc
tic
%% Parameters, etc.
Temp = linspace(350+273,550+273,10);
PN2 = 200/4*0.8;
PH2 = 200*3/4*0.8;
PNH3 = 0.2*200;

norm = 10^6; % just to keep values around unity
% n, nd, h, hd, r1, r1_rev, r2, r2_rev, r3, r3_rev, dd, d
Arow = 1.46*10^13*ones(8,1);
Act = [2.39*10^4; 1.46*10^14; 3.89*10^6; Arow; 2.46*10^3; 1.46*10^13; 1.46*10^13; 1.46*10^13; 1.46*10^13; 1.46*10^13]; % 1/s
Ea = [0.4; 1.2; 0; 1.1; 1.1; 1.1; 1.3; 1.46; 1.2; 0.75; 0.85; 0; 0.94; 0.68; 0.41; 0.15; 0.21]; %eV
kB = 8.617333262*10^-5; % eV/K

% % turn off diffusion
% Act(end-4:end) = 0;

tspan = [0 10^20];
delta = 10^-3; % differential change for calculation of reaction sensitivities
eps = 10^-2; % all adsorbates start at 1% coverage
%% Output arrays
S = zeros(length(Temp),length(Act))';

A = zeros(1,length(Temp)); B = A; C = A; D = A; H = A; O = A;
AA = A; AO = A; AB = A; AC = A; AD = A; AH = A;
BH = A; CH = A; DH = A; HO = A; HH = A;
BB = A; BC = A; BD = A; BO = A;
CC = A; CD = A; CO = A;
DD = A; DO = A;
OO = A;

r = A;

ZN = A; ZH = A; ZR1 = A; ZR2 = A; ZR3 = A; ZM = A;
%% Solvers and Nested Loops
for w = 1:length(Temp)
    disp(w)
    toc
    T = Temp(w);
    k00 = Act.*exp(-Ea/kB./T)/norm;
    k0 = k00;
    k0(1) = k00(1)*PN2;
    k0(3) = k00(3)*PH2;
    k0(12) = k00(12)*PNH3;
    for m = 1:1
        for j = 1:length(k0)
            K = zeros(2,length(k0));
            K(1,:) = k0; K(2,:) = k0;
            K(2,j) = k0(j)*(1+delta);
            for p = 1:2
                %% Rate Constant assignments
                % n, nd, h, hd, r1, r1_rev, r2, r2_rev, r3, r3_rev, dd, d
                kn = K(p,1);
                knd = K(p,2);
                kh = K(p,3);
                khd = K(p,4);
                kr1 = K(p,5);
                kr1_rev = K(p,6);
                kr2 = K(p,7);
                kr2_rev = K(p,8);
                kr3 = K(p,9);
                kr3_rev = K(p,10);
                kdd = K(p,11);
                kd = K(p,12);
                k_diff_n = K(p,13);
                k_diff_nh = K(p,14);
                k_diff_nh2 = K(p,15);
                k_diff_nh3 = K(p,16);
                k_diff_h = K(p,17);
                %% Sites
                % a = N*, b = NH*, c = NH2*, d = NH3*, h = H*
                d = @(x) x(1); c = @(x) x(2); b = @(x) x(3); a = @(x) x(4); h = @(x) x(5);
                dd = @(x) x(6); cd = @(x) x(7); bd = @(x) x(8); ad = @(x) x(9); hd = @(x) x(10);
                cc = @(x) x(11); bc = @(x) x(12); ac = @(x) x(13); hc = @(x) x(14);
                bb = @(x) x(15); ab = @(x) x(16); hb = @(x) x(17);
                aa = @(x) x(18); ha = @(x) x(19);
                hh = @(x) x(20);
                
                o = @(x) 1 - a(x) - b(x) - c(x) - d(x) - h(x); % overall site balance
                do = @(x) d(x) - dd(x) - cd(x) - bd(x) - ad(x) - hd(x); % d site balance
                co = @(x) c(x) - cd(x) - cc(x) - bc(x) - ac(x) - hc(x); % c site balance
                bo = @(x) b(x) - bd(x) - bc(x) - bb(x) - ab(x) - hb(x); % b site balance
                ao = @(x) a(x) - ad(x) - ac(x) - ab(x) - aa(x) - ha(x); % a site balance
                ho = @(x) h(x) - hd(x) - hc(x) - hb(x) - ha(x) - hh(x); % h site balance
                oo = @(x) o(x) - do(x) - co(x) - bo(x) - ao(x) - ho(x); % o site balance    
                %% N2 adsorption
                % d, c, b, a, h
                n_d = @(x) 0; n_c = @(x) 0; n_b = @(x) 0; n_a = @(x) 4*kn*oo(x); n_h = @(x) 0;
                n_dd = @(x) 0; n_cd = @(x) 0; n_bd = @(x) 0; n_ad = @(x) 3*kn.*oo(x).*do(x)./o(x); n_hd = @(x) 0;
                n_cc = @(x) 0; n_bc = @(x) 0; n_ac = @(x) 3*kn.*oo(x).*co(x)./o(x); n_hc = @(x) 0;
                n_bb = @(x) 0; n_ab = @(x) 3*kn.*oo(x).*bo(x)./o(x); n_hb = @(x) 0;
                n_aa = @(x) kn*oo(x).*(1 + 6*ao(x)./o(x)); n_ha = @(x) 3*kn.*oo(x).*ho(x)./o(x);
                n_hh = @(x) 0;
                %% N2 desorption
                nd_d = @(x) 0; nd_c = @(x) 0; nd_b = @(x) 0; nd_a = @(x) -4*knd*aa(x); nd_h = @(x) 0;
                nd_dd = @(x) 0; nd_cd = @(x) 0; nd_bd = @(x) 0; nd_ad = @(x) -3*knd.*aa(x).*ad(x)./a(x); nd_hd = @(x) 0;
                nd_cc = @(x) 0; nd_bc = @(x) 0; nd_ac = @(x) -3*knd.*aa(x).*ac(x)./a(x); nd_hc = @(x) 0;
                nd_bb = @(x) 0; nd_ab = @(x) -3*knd.*aa(x).*ab(x)./a(x); nd_hb = @(x) 0;
                nd_aa = @(x) -knd*aa(x).*(1 + 6*aa(x)./a(x)); nd_ha = @(x) -3*knd.*aa(x).*ha(x)./a(x);
                nd_hh = @(x) 0;
                %% H2 adsorption
                h_d = @(x) 0; h_c = @(x) 0; h_b = @(x) 0; h_a = @(x) 0; h_h = @(x) 4*kh*oo(x);
                h_dd = @(x) 0; h_cd = @(x) 0; h_bd = @(x) 0; h_ad = @(x) 0; h_hd = @(x) 3*kh*oo(x).*do(x)./o(x);
                h_cc = @(x) 0; h_bc = @(x) 0; h_ac = @(x) 0; h_hc = @(x) 3*kh*oo(x).*co(x)./o(x);
                h_bb = @(x) 0; h_ab = @(x) 0; h_hb = @(x) 3*kh*oo(x).*bo(x)./o(x);
                h_aa = @(x) 0; h_ha = @(x) 3*kh.*oo(x).*ao(x)./o(x);
                h_hh = @(x) kh*oo(x).*(1+6.*ho(x)./o(x));
                %% H2 desorption
                hd_d = @(x) 0; hd_c = @(x) 0; hd_b = @(x) 0; hd_a = @(x) 0; hd_h = @(x) -4*khd*hh(x);
                hd_dd = @(x) 0; hd_cd = @(x) 0; hd_bd = @(x) 0; hd_ad = @(x) 0; hd_hd = @(x) -3*khd*hh(x).*hd(x)./h(x);
                hd_cc = @(x) 0; hd_bc = @(x) 0; hd_ac = @(x) 0; hd_hc = @(x) -3*khd*hh(x).*hc(x)./h(x);
                hd_bb = @(x) 0; hd_ab = @(x) 0; hd_hb = @(x) -3*khd*hh(x).*hb(x)./h(x);
                hd_aa = @(x) 0; hd_ha = @(x) -3*khd.*hh(x).*ha(x)./h(x);
                hd_hh = @(x) -khd*hh(x).*(1+6.*hh(x)./h(x));
                %% N+H (ha -> ob)
                r1_d = @(x) 0; r1_c = @(x) 0; r1_b = @(x) 4*kr1*ha(x); r1_a = @(x) -4*kr1*ha(x); r1_h = @(x) -4*kr1*ha(x);
                r1_dd = @(x) 0; r1_cd = @(x) 0; r1_bd = @(x) 3*kr1*ha(x).*ad(x)./a(x); r1_ad = @(x) -3*kr1*ha(x).*ad(x)./a(x); r1_hd = @(x) -3*kr1*ha(x).*hd(x)./h(x);
                r1_cc = @(x) 0; r1_bc = @(x) 3*kr1*ha(x).*ac(x)./a(x); r1_ac = @(x) -3*kr1*ha(x).*ac(x)./a(x); r1_hc = @(x) -3*kr1*ha(x).*hc(x)./h(x);
                r1_bb = @(x) 6*kr1.*ha(x).*ab(x)./a(x); r1_ab = @(x) 3*kr1*ha(x).*(aa(x)./a(x) - ab(x)./a(x)); r1_hb = @(x) 3*kr1*ha(x).*(ha(x)./a(x) - hb(x)./h(x));
                r1_aa = @(x) -6*kr1.*ha(x).*aa(x)./a(x); r1_ha = @(x) -kr1.*ha(x).*(1+3.*ha(x)./h(x)+3.*ha(x)./a(x));
                r1_hh = @(x) -6*kr1.*ha(x).*hh(x)./h(x);
                %% NH decomp. (ob -> ha)
                r1_rev_d = @(x) 0; r1_rev_c = @(x) 0; r1_rev_b = @(x) -4*kr1_rev*bo(x); r1_rev_a = @(x) 4*kr1_rev*bo(x); r1_rev_h = @(x) 4*kr1_rev*bo(x);
                r1_rev_dd = @(x) 0; r1_rev_cd = @(x) 0; r1_rev_bd = @(x) -3*kr1_rev*bo(x).*bd(x)./b(x); r1_rev_ad = @(x) 3*kr1_rev*bo(x).*bd(x)./b(x); r1_rev_hd = @(x) 3*kr1_rev*bo(x).*do(x)./o(x);
                r1_rev_cc = @(x) 0; r1_rev_bc = @(x) -3*kr1_rev*bo(x).*bc(x)./b(x); r1_rev_ac = @(x) 3*kr1_rev*bo(x).*bc(x)./b(x); r1_rev_hc = @(x) 3*kr1_rev*bo(x).*co(x)./o(x);
                r1_rev_bb = @(x) -6*kr1_rev.*bo(x).*bb(x)./b(x); r1_rev_ab = @(x) 3*kr1_rev*bo(x).*(bb(x)./b(x) - ab(x)./b(x)); r1_rev_hb = @(x) 3*kr1_rev*bo(x).*(bo(x)./o(x) - hb(x)./b(x));
                r1_rev_aa = @(x) 6*kr1_rev.*bo(x).*ab(x)./b(x); r1_rev_ha = @(x) kr1_rev.*bo(x).*(1+3.*hb(x)./b(x)+3.*ao(x)./o(x));
                r1_rev_hh = @(x) 6*kr1_rev.*bo(x).*ho(x)./o(x);  
                %% NH+H (hb -> oc)
                r2_d = @(x) 0; r2_c = @(x) 4*kr2*hb(x); r2_b = @(x) -4*kr2*hb(x); r2_a = @(x) 0; r2_h = @(x) -4*kr2*hb(x);
                r2_dd = @(x) 0; r2_cd = @(x) 3*kr2*hb(x).*bd(x)./b(x); r2_bd = @(x) -3*kr2*hb(x).*bd(x)./b(x); r2_ad = @(x) 0; r2_hd = @(x) -3*kr2*hb(x).*hd(x)./h(x);
                r2_cc = @(x) 6*kr2.*hb(x).*bc(x)./b(x); r2_bc = @(x) 3*kr2.*hb(x).*(bb(x)./b(x) - bc(x)./b(x)); r2_ac = @(x) 3*kr2.*hb(x).*ab(x)./b(x); r2_hc = @(x) 3*kr2.*hb(x).*(hb(x)./b(x) - hc(x)./h(x));
                r2_bb = @(x) -6*kr2.*hb(x).*bb(x)./b(x); r2_ab = @(x) -3*kr2*hb(x).*ab(x)./b(x); r2_hb = @(x) -kr2.*hb(x).*(1 + 3.*hb(x)./h(x) + 3*hb(x)./b(x));
                r2_aa = @(x) 0; r2_ha = @(x) -3*kr2.*hb(x).*ha(x)./h(x);
                r2_hh = @(x) -6*kr2.*hb(x).*hh(x)./h(x);
                %% NH2 -> NH+H (oc -> hb)
                r2_rev_d = @(x) 0; r2_rev_c = @(x) -4*kr2_rev*co(x); r2_rev_b = @(x) 4*kr2_rev*co(x); r2_rev_a = @(x) 0; r2_rev_h = @(x) 4*kr2_rev*co(x);
                r2_rev_dd = @(x) 0; r2_rev_cd = @(x) -3*kr2_rev*co(x).*cd(x)./c(x); r2_rev_bd = @(x) 3*kr2_rev*co(x).*cd(x)./c(x); r2_rev_ad = @(x) 0; r2_rev_hd = @(x) 3*kr2_rev*co(x).*do(x)./o(x);
                r2_rev_cc = @(x) -6*kr2_rev.*co(x).*cc(x)./c(x); r2_rev_bc = @(x) 3*kr2_rev.*co(x).*(cc(x)./c(x) - bc(x)./c(x)); r2_rev_ac = @(x) -3*kr2_rev.*co(x).*ac(x)./c(x); r2_rev_hc = @(x) 3*kr2_rev.*co(x).*(co(x)./o(x) - hc(x)./c(x));
                r2_rev_bb = @(x) 6*kr2_rev.*co(x).*bc(x)./c(x); r2_rev_ab = @(x) 3*kr2_rev*co(x).*ac(x)./c(x); r2_rev_hb = @(x) kr2_rev.*co(x).*(1 + 3.*co(x)./o(x) + 3*co(x)./c(x));
                r2_rev_aa = @(x) 0; r2_rev_ha = @(x) 3*kr2_rev.*co(x).*ao(x)./o(x);
                r2_rev_hh = @(x) 6*kr2_rev.*co(x).*ho(x)./o(x);
                %% NH2 + H (ch -> do)
                r3_d = @(x) 4*kr3*hc(x); r3_c = @(x) -4*kr3*hc(x); r3_b = @(x) 0; r3_a = @(x) 0; r3_h = @(x) -4*kr3*hc(x);
                r3_dd = @(x) 6*kr3.*hc(x).*cd(x)./c(x); r3_cd = @(x) 3*kr3*hc(x).*(cc(x)./c(x) - cd(x)./c(x)); r3_bd = @(x) 3*kr3*hc(x).*bc(x)./c(x); r3_ad = @(x) 3*kr3*hc(x).*ac(x)./c(x); r3_hd = @(x) 3*kr3*hc(x).*(hc(x)./c(x) - hd(x)./h(x));
                r3_cc = @(x) -6*kr3.*hc(x).*cc(x)./c(x); r3_bc = @(x) -3*kr3*hc(x).*bc(x)./c(x); r3_ac = @(x) -3*kr3*hc(x).*ac(x)./c(x); r3_hc = @(x) -kr3.*hc(x).*(1+3*hc(x)./h(x)+3*hc(x)./c(x));
                r3_bb = @(x) 0; r3_ab = @(x) 0; r3_hb = @(x) -3*kr3*hc(x).*hb(x)./h(x);
                r3_aa = @(x) 0; r3_ha = @(x) -3*kr3.*hc(x).*ha(x)./h(x);
                r3_hh = @(x) -6*kr3.*hc(x).*hh(x)./h(x);              
                %% NH3 decomp. (do -> ch)
                r3_rev_d = @(x) -4*kr3_rev*do(x); r3_rev_c = @(x) 4*kr3_rev*do(x); r3_rev_b = @(x) 0; r3_rev_a = @(x) 0; r3_rev_h = @(x) 4*kr3_rev*do(x);
                r3_rev_dd = @(x) -6*kr3_rev.*do(x).*dd(x)./d(x); r3_rev_cd = @(x) 3*kr3_rev*do(x).*(dd(x)./d(x) - cd(x)./d(x)); r3_rev_bd = @(x) -3*kr3_rev*do(x).*bd(x)./d(x); r3_rev_ad = @(x) -3*kr3_rev*do(x).*ad(x)./d(x); r3_rev_hd = @(x) 3*kr3_rev*do(x).*(do(x)./o(x) - hd(x)./d(x));
                r3_rev_cc = @(x) 6*kr3_rev.*do(x).*cd(x)./d(x); r3_rev_bc = @(x) 3*kr3_rev*do(x).*bd(x)./d(x); r3_rev_ac = @(x) 3*kr3_rev*do(x).*ad(x)./d(x); r3_rev_hc = @(x) kr3_rev.*do(x).*(1+3*hd(x)./d(x)+3*co(x)./o(x));
                r3_rev_bb = @(x) 0; r3_rev_ab = @(x) 0; r3_rev_hb = @(x) 3*kr3_rev*do(x).*bo(x)./o(x);
                r3_rev_aa = @(x) 0; r3_rev_ha = @(x) 3*kr3_rev.*do(x).*ao(x)./o(x);
                r3_rev_hh = @(x) 6*kr3_rev.*do(x).*ho(x)./o(x);               
                %% NH3 desorp. (d -> o)
                dd_d = @(x) -kdd.*d(x); dd_c = @(x) 0; dd_b = @(x) 0; dd_a = @(x) 0; dd_h = @(x) 0;
                dd_dd = @(x) -2*kdd.*dd(x); dd_cd = @(x) -kdd.*cd(x); dd_bd = @(x) -kdd.*bd(x); dd_ad = @(x) -kdd.*ad(x); dd_hd = @(x) -kdd.*hd(x);
                dd_cc = @(x) 0; dd_bc = @(x) 0; dd_ac = @(x) 0; dd_hc = @(x) 0;
                dd_bb = @(x) 0; dd_ab = @(x) 0; dd_hb = @(x) 0;
                dd_aa = @(x) 0; dd_ha = @(x) 0;
                dd_hh = @(x) 0;
                %% NH3 adsorp. (o -> d)
                d_d = @(x) kd.*o(x); d_c = @(x) 0; d_b = @(x) 0; d_a = @(x) 0; d_h = @(x) 0;
                d_dd = @(x) 2*kd.*do(x); d_cd = @(x) kd.*co(x); d_bd = @(x) kd.*bo(x); d_ad = @(x) kd.*ao(x); d_hd = @(x) kd.*ho(x);
                d_cc = @(x) 0; d_bc = @(x) 0; d_ac = @(x) 0; d_hc = @(x) 0;
                d_bb = @(x) 0; d_ab = @(x) 0; d_hb = @(x) 0;
                d_aa = @(x) 0; d_ha = @(x) 0;
                d_hh = @(x) 0;
                %% N diffusion (ao -> oa)
                diff_d = @(x) 0; diff_c = @(x) 0; diff_b = @(x) 0; diff_a = @(x) 0; diff_h = @(x) 0;
                diff_dd = @(x) 0; diff_cd = @(x) 0; diff_bd = @(x) 0; diff_ad = @(x) 3*k_diff_n*ao(x).*(do(x)./o(x) - ad(x)./a(x)); diff_hd = @(x) 0;
                diff_cc = @(x) 0; diff_bc = @(x) 0; diff_ac = @(x) 3*k_diff_n*ao(x).*(co(x)./o(x) - ac(x)./a(x)); diff_hc = @(x) 0;
                diff_bb = @(x) 0; diff_ab = @(x) 3*k_diff_n*ao(x).*(bo(x)./o(x) - ab(x)./a(x)); diff_hb = @(x) 0;
                diff_aa = @(x) 6*k_diff_n*ao(x).*(ao(x)./o(x) - aa(x)./a(x)); diff_ha = @(x) 3*k_diff_n*ao(x).*(ho(x)./o(x) - ha(x)./a(x));
                diff_hh = @(x) 0;
                %% NH diffusion (bo -> ob)
                diff2_d = @(x) 0; diff2_c = @(x) 0; diff2_b = @(x) 0; diff2_a = @(x) 0; diff2_h = @(x) 0;
                diff2_dd = @(x) 0; diff2_cd = @(x) 0; diff2_bd = @(x) 3*k_diff_nh*bo(x).*(do(x)./o(x) - bd(x)./b(x)); diff2_ad = @(x) 0; diff2_hd = @(x) 0;
                diff2_cc = @(x) 0; diff2_bc = @(x) 3*k_diff_nh*bo(x).*(co(x)./o(x) - bc(x)./b(x)); diff2_ac = @(x) 0; diff2_hc = @(x) 0;
                diff2_bb = @(x) 6*k_diff_nh*bo(x).*(bo(x)./o(x) - bb(x)./b(x)); diff2_ab = @(x) 3*k_diff_nh*bo(x).*(ao(x)./o(x) - ab(x)./b(x)); diff2_hb = @(x) 0;
                diff2_aa = @(x) 0; diff2_ha = @(x) 0;
                diff2_hh = @(x) 0;
                %% NH2 diffusion (co -> oc)
                diff3_d = @(x) 0; diff3_c = @(x) 0; diff3_b = @(x) 0; diff3_a = @(x) 0; diff3_h = @(x) 0;
                diff3_dd = @(x) 0; diff3_cd = @(x) 3*k_diff_nh2*co(x).*(do(x)./o(x) - cd(x)./c(x)); diff3_bd = @(x) 0; diff3_ad = @(x) 0; diff3_hd = @(x) 0;
                diff3_cc = @(x) 6*k_diff_nh2*co(x).*(co(x)./o(x) - cc(x)./c(x)); diff3_bc = @(x) 3*k_diff_nh2*co(x).*(bo(x)./o(x) - bc(x)./c(x)); diff3_ac = @(x) 3*k_diff_nh2*co(x).*(ao(x)./o(x) - ac(x)./c(x)); diff3_hc = @(x) 3*k_diff_nh2*co(x).*(ho(x)./o(x) - hc(x)./c(x));
                diff3_bb = @(x) 0; diff3_ab = @(x) 0; diff3_hb = @(x) 0;
                diff3_aa = @(x) 0; diff3_ha = @(x) 0;
                diff3_hh = @(x) 0;
                %% NH3 diffusion (do -> od)
                diff4_d = @(x) 0; diff4_c = @(x) 0; diff4_b = @(x) 0; diff4_a = @(x) 0; diff4_h = @(x) 0;
                diff4_dd = @(x) 6*k_diff_nh3*do(x).*(do(x)./o(x) - dd(x)./d(x)); diff4_cd = @(x) 3*k_diff_nh3*do(x).*(co(x)./o(x) - cd(x)./d(x)); diff4_bd = @(x) 3*k_diff_nh3*do(x).*(bo(x)./o(x) - bd(x)./d(x)); diff4_ad = @(x) 3*k_diff_nh3*do(x).*(ao(x)./o(x) - ad(x)./d(x)); diff4_hd = @(x) 3*k_diff_nh3*do(x).*(ho(x)./o(x) - hd(x)./d(x));
                diff4_cc = @(x) 0; diff4_bc = @(x) 0; diff4_ac = @(x) 0; diff4_hc = @(x) 0;
                diff4_bb = @(x) 0; diff4_ab = @(x) 0; diff4_hb = @(x) 0;
                diff4_aa = @(x) 0; diff4_ha = @(x) 0;
                diff4_hh = @(x) 0;
                %% H diffusion (ho -> oh)
                diff5_d = @(x) 0; diff5_c = @(x) 0; diff5_b = @(x) 0; diff5_a = @(x) 0; diff5_h = @(x) 0;
                diff5_dd = @(x) 0; diff5_cd = @(x) 0; diff5_bd = @(x) 0; diff5_ad = @(x) 0; diff5_hd = @(x) 3*k_diff_h*ho(x).*(do(x)./o(x) - hd(x)./h(x));
                diff5_cc = @(x) 0; diff5_bc = @(x) 0; diff5_ac = @(x) 0; diff5_hc = @(x) 3*k_diff_h*ho(x).*(co(x)./o(x) - hc(x)./h(x));
                diff5_bb = @(x) 0; diff5_ab = @(x) 0; diff5_hb = @(x) 3*k_diff_h*ho(x).*(bo(x)./o(x) - hb(x)./h(x));
                diff5_aa = @(x) 0; diff5_ha = @(x) 3*k_diff_h*ho(x).*(ao(x)./o(x) - ha(x)./h(x));
                diff5_hh = @(x) 6*k_diff_h*ho(x).*(ho(x)./o(x) - hh(x)./h(x));
                %% ODE Matrix 
                F = @(t,x) [n_d(x) + nd_d(x) + h_d(x) + hd_d(x) + r1_d(x) + r1_rev_d(x) + r2_d(x) + r2_rev_d(x) + r3_d(x) + r3_rev_d(x) + dd_d(x) + d_d(x) + diff_d(x) + diff2_d(x) + diff3_d(x) + diff4_d(x) + diff5_d(x);
                    n_c(x) + nd_c(x) + h_c(x) + hd_c(x) + r1_c(x) + r1_rev_c(x) + r2_c(x) + r2_rev_c(x) + r3_c(x) + r3_rev_c(x) + dd_c(x) + d_c(x) + diff_c(x) + diff2_c(x) + diff3_c(x) + diff4_c(x) + diff5_c(x);
                    n_b(x) + nd_b(x) + h_b(x) + hd_b(x) + r1_b(x) + r1_rev_b(x) + r2_b(x) + r2_rev_b(x) + r3_b(x) + r3_rev_b(x) + dd_b(x) + d_b(x) + diff_b(x) + diff2_b(x) + diff3_b(x) + diff4_b(x) + diff5_b(x);
                    n_a(x) + nd_a(x) + h_a(x) + hd_a(x) + r1_a(x) + r1_rev_a(x) + r2_a(x) + r2_rev_a(x) + r3_a(x) + r3_rev_a(x) + dd_a(x) + d_a(x) + diff_a(x) + diff2_a(x) + diff3_a(x) + diff4_a(x) + diff5_a(x);
                    n_h(x) + nd_h(x) + h_h(x) + hd_h(x) + r1_h(x) + r1_rev_h(x) + r2_h(x) + r2_rev_h(x) + r3_h(x) + r3_rev_h(x) + dd_h(x) + d_h(x)  + diff_h(x) + diff2_h(x) + diff3_h(x) + diff4_h(x) + diff5_h(x);
                    n_dd(x) + nd_dd(x) + h_dd(x) + hd_dd(x) + r1_dd(x) + r1_rev_dd(x) + r2_dd(x) + r2_rev_dd(x) + r3_dd(x) + r3_rev_dd(x) + dd_dd(x) + d_dd(x) + diff_dd(x) + diff2_dd(x) + diff3_dd(x) + diff4_dd(x) + diff5_dd(x);
                    n_cd(x) + nd_cd(x) + h_cd(x) + hd_cd(x) + r1_cd(x) + r1_rev_cd(x) + r2_cd(x) + r2_rev_cd(x) + r3_cd(x) + r3_rev_cd(x) + dd_cd(x) + d_cd(x) + diff_cd(x) + diff2_cd(x) + diff3_cd(x) + diff4_cd(x) + diff5_cd(x);
                    n_bd(x) + nd_bd(x) + h_bd(x) + hd_bd(x) + r1_bd(x) + r1_rev_bd(x) + r2_bd(x) + r2_rev_bd(x) + r3_bd(x) + r3_rev_bd(x) + dd_bd(x) + d_bd(x) + diff_bd(x) + diff2_bd(x) + diff3_bd(x) + diff4_bd(x) + diff5_bd(x);
                    n_ad(x) + nd_ad(x) + h_ad(x) + hd_ad(x) + r1_ad(x) + r1_rev_ad(x) + r2_ad(x) + r2_rev_ad(x) + r3_ad(x) + r3_rev_ad(x) + dd_ad(x) + d_ad(x) + diff_ad(x) + diff2_ad(x) + diff3_ad(x) + diff4_ad(x) + diff5_ad(x);
                    n_hd(x) + nd_hd(x) + h_hd(x) + hd_hd(x) + r1_hd(x) + r1_rev_hd(x) + r2_hd(x) + r2_rev_hd(x) + r3_hd(x) + r3_rev_hd(x) + dd_hd(x) + d_hd(x) + diff_hd(x) + diff2_hd(x) + diff3_hd(x) + diff4_hd(x) + diff5_hd(x);
                    n_cc(x) + nd_cc(x) + h_cc(x) + hd_cc(x) + r1_cc(x) + r1_rev_cc(x) + r2_cc(x) + r2_rev_cc(x) + r3_cc(x) + r3_rev_cc(x) + dd_cc(x) + d_cc(x) + diff_cc(x) + diff2_cc(x) + diff3_cc(x) + diff4_cc(x) + diff5_cc(x);
                    n_bc(x) + nd_bc(x) + h_bc(x) + hd_bc(x) + r1_bc(x) + r1_rev_bc(x) + r2_bc(x) + r2_rev_bc(x) + r3_bc(x) + r3_rev_bc(x) + dd_bc(x) + d_bc(x) + diff_bc(x) + diff2_bc(x) + diff3_bc(x) + diff4_bc(x) + diff5_bc(x);
                    n_ac(x) + nd_ac(x) + h_ac(x) + hd_ac(x) + r1_ac(x) + r1_rev_ac(x) + r2_ac(x) + r2_rev_ac(x) + r3_ac(x) + r3_rev_ac(x) + dd_ac(x) + d_ac(x) + diff_ac(x) + diff2_ac(x) + diff3_ac(x) + diff4_ac(x) + diff5_ac(x);
                    n_hc(x) + nd_hc(x) + h_hc(x) + hd_hc(x) + r1_hc(x) + r1_rev_hc(x) + r2_hc(x) + r2_rev_hc(x) + r3_hc(x) + r3_rev_hc(x) + dd_hc(x) + d_hc(x) + diff_hc(x) + diff2_hc(x) + diff3_hc(x) + diff4_hc(x) + diff5_hc(x);
                    n_bb(x) + nd_bb(x) + h_bb(x) + hd_bb(x) + r1_bb(x) + r1_rev_bb(x) + r2_bb(x) + r2_rev_bb(x) + r3_bb(x) + r3_rev_bb(x) + dd_bb(x) + d_bb(x) + diff_bb(x) + diff2_bb(x) + diff3_bb(x) + diff4_bb(x) + diff5_bb(x);
                    n_ab(x) + nd_ab(x) + h_ab(x) + hd_ab(x) + r1_ab(x) + r1_rev_ab(x) + r2_ab(x) + r2_rev_ab(x) + r3_ab(x) + r3_rev_ab(x) + dd_ab(x) + d_ab(x) + diff_ab(x) + diff2_ab(x) + diff3_ab(x) + diff4_ab(x) + diff5_ab(x);
                    n_hb(x) + nd_hb(x) + h_hb(x) + hd_hb(x) + r1_hb(x) + r1_rev_hb(x) + r2_hb(x) + r2_rev_hb(x) + r3_hb(x) + r3_rev_hb(x) + dd_hb(x) + d_hb(x) + diff_hb(x) + diff2_hb(x) + diff3_hb(x) + diff4_hb(x) + diff5_hb(x);
                    n_aa(x) + nd_aa(x) + h_aa(x) + hd_aa(x) + r1_aa(x) + r1_rev_aa(x) + r2_aa(x) + r2_rev_aa(x) + r3_aa(x) + r3_rev_aa(x) + dd_aa(x) + d_aa(x) + diff_aa(x) + diff2_aa(x) + diff3_aa(x) + diff4_aa(x) + diff5_aa(x);
                    n_ha(x) + nd_ha(x) + h_ha(x) + hd_ha(x) + r1_ha(x) + r1_rev_ha(x) + r2_ha(x) + r2_rev_ha(x) + r3_ha(x) + r3_rev_ha(x) + dd_ha(x) + d_ha(x) + diff_ha(x) + diff2_ha(x) + diff3_ha(x) + diff4_ha(x) + diff5_ha(x);
                    n_hh(x) + nd_hh(x) + h_hh(x) + hd_hh(x) + r1_hh(x) + r1_rev_hh(x) + r2_hh(x) + r2_rev_hh(x) + r3_hh(x) + r3_rev_hh(x) + dd_hh(x) + d_hh(x) + diff_hh(x) + diff2_hh(x) + diff3_hh(x) + diff4_hh(x) + diff5_hh(x)];                         
                C00 = ones(1,5)*eps;
                C01 = ones(1,15)*eps^2;
                C0 = [C00 C01];
                %% ODE15s
                [T,Y] = ode23s(F,tspan,C0); t = T;
                d = Y(:,1); c = Y(:,2); b = Y(:,3); a = Y(:,4); h = Y(:,5);
                dd = Y(:,6); cd = Y(:,7); bd = Y(:,8); ad = Y(:,9); hd = Y(:,10);
                cc = Y(:,11); bc = Y(:,12); ac = Y(:,13); hc = Y(:,14);
                bb = Y(:,15); ab = Y(:,16); hb = Y(:,17);
                aa = Y(:,18); ha = Y(:,19);
                hh = Y(:,20);
                
                o = 1 - a - b - c - d - h; % overall site balance
                do = d - dd - cd - bd - ad - hd; % d site balance
                co = c - cd - cc - bc - ac - hc; % c site balance
                bo = b - bd - bc - bb - ab - hb; % b site balance
                ao = a - ad - ac - ab - aa - ha; % a site balance
                ho = h - hd - hc - hb - ha - hh; % h site balance
                oo = o - do - co - bo - ao - ho; % o site balance
                
                em = round(length(oo)*0.10);
                
                % Elementary step reversibilities
                zN = (knd.*mean(aa(end-em:end)))./(kn.*mean(oo(end-em:end)));
                zH = (khd.*mean(hh(end-em:end)))./(kh.*mean(oo(end-em:end)));
                zr1 = (kr1_rev.*mean(bo(end-em:end)))./(kr1.*mean(ha(end-em:end)));
                zr2 = (kr2_rev.*mean(co(end-em:end)))./(kr2.*mean(hb(end-em:end)));
                zr3 = (kr3_rev.*mean(do(end-em:end)))./(kr3.*mean(hc(end-em:end)));
                zM = (kd.*mean(o(end-em:end)))./(kdd.*mean(d(end-em:end)));
                
                R(p,j) = kn.*mean(oo(end-em:end)) - knd.*mean(aa(end-em:end));
               
%                 if mod(w,round(length(P))/10) == 0 && j == length(k0) && m == 1 && p == 1
%                     figure
%                     set(gca,'fontweight','bold','fontsize',11,'box','on');
%                     hold on
%                     set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
%                     title(num2str(w))
%                     plot(log10(t),a,'linewidth',1)
%                     plot(log10(t),b,'linewidth',1)
%                     plot(log10(t),c,'linewidth',1)
%                     plot(log10(t),d,'linewidth',1)
%                     plot(log10(t),h,'linewidth',1)
%                     plot(log10(t),o,'linewidth',1)
%                     legend('a','b','c','d','h','o')
%                 end
            end
            %% Extract coverages, DORC
            S(j,w) = (R(2,j) - R(1,j))./(delta*R(1,j));
            r(m,w) = kdd.*mean(d(end-round(length(d)*0.10):end)) - kd.*mean(o(end-round(length(o)*0.10):end));
            O(m,w) = mean(o(end-em:end));A(m,w) = mean(a(end-em:end));B(m,w) = mean(b(end-em:end));C(m,w) = mean(c(end-em:end));D(m,w) = mean(d(end-em:end));H(m,w) = mean(h(end-em:end));
            AH(m,w) = mean(ha(end-em:end));BH(m,w) = mean(hb(end-em:end));CH(m,w) = mean(hc(end-em:end));DH(m,w) = mean(hd(end-em:end));HH(m,w) = mean(hh(end-em:end)); HO(m,w) = mean(ho(end-em:end));
            OO(m,w) = mean(oo(end-em:end));AO(m,w) = mean(ao(end-em:end));BO(m,w) = mean(bo(end-em:end));CO(m,w) = mean(co(end-em:end));DO(m,w) = mean(do(end-em:end));
            AA(m,w) = mean(aa(end-em:end));AB(m,w) = mean(ab(end-em:end));AC(m,w) = mean(ac(end-em:end));AD(m,w) = mean(ad(end-em:end));
            BB(m,w) = mean(bb(end-em:end));BC(m,w) = mean(bc(end-em:end));BD(m,w) = mean(bd(end-em:end));
            CC(m,w) = mean(cc(end-em:end));CD(m,w) = mean(cd(end-em:end));
            DD(m,w) = mean(dd(end-em:end));
            
            ZN(w) = zN; ZH(w) = zH; ZR1(w) = zr1; ZR2(w) = zr2; ZR3(w) = zr3; ZM(w) = zM;
        end
    end
end
%% Extract results
mu_AA = AA./A./A; mu_AB = AB./A./B; mu_AC = AC./A./C; mu_AD = AD./A./D; mu_AO = AO./A./O; mu_AH = AH./A./H;
mu_BB = BB./B./B; mu_BC = BC./B./C; mu_BD = BD./B./D; mu_BO = BO./B./O; mu_BH = BH./B./H;
mu_CC = CC./C./C; mu_CD = CD./C./D; mu_CO = CO./C./O; mu_CH = CH./C./H;
mu_DD = DD./D./D; mu_DO = DO./D./O; mu_DH = DH./D./H;
mu_HH = HH./H./H; mu_HO = HO./H./O;
mu_OO = OO./O./O; 

App = (Ea)'*S*1.60218e-19*6.0221409e+23/1000; % Apparent activation energy, change units from eV/atom to kJ/mol
%% DORC sum to unity
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(Temp,sum(S,1),'k-','linewidth',1)
ylim([0 2]);
title('DORC check');
%% Apparent activation energy
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(Temp,App,'k-','linewidth',1)
xlabel('T / K');
ylabel('E_{A,app} / kJ mol^{-1}');
%% Apparent activation energy -- elementary step contributions
% N2 -sorption, H2 desorption, N* hydrogenation are the only contributors to E_A,app
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(Temp,S.*Ea*1.60218e-19*6.0221409e+23/1000,'linewidth',1)
xlabel('T / K');
ylabel('E_{A,app} / kJ mol^{-1}');
legend('1','-1','6','-6','2','-2','3','-3','4','-4','5','-5','diff');
legend boxoff
%% DORC's
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(Temp,S(1,:)+S(2,:),'linewidth',1.5)
plot(Temp,S(3,:)+S(4,:),'linewidth',1.5)
plot(Temp,S(5,:)+S(6,:),'linewidth',1.5)
plot(Temp,S(7,:)+S(8,:),'linewidth',1.5)
plot(Temp,S(9,:)+S(10,:),'linewidth',1.5)
plot(Temp,S(11,:)+S(12,:),'linewidth',1.5)
plot(Temp,S(13,:),'k-.','linewidth',1.5)
xlabel('T / K');
legend('X_{RC,N_{2} ads}','X_{RC,H_{2} ads}','X_{RC,r1}','X_{RC,r2}','X_{RC,r3}','X_{RC,NH_{3} des}','diff');
legend boxoff
legend('location','northwest');
ylim([-0.2 1.2])
%% Coverages
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(Temp,log10(A),'linewidth',1.5)
plot(Temp,log10(H),'linewidth',1.5)
plot(Temp,log10(D),'linewidth',1.5)
plot(Temp,log10(O),'linewidth',1.5)
plot(Temp,log10(B),'linewidth',1.5)
plot(Temp,log10(C),'linewidth',1.5)
xlabel('T / K');
legend('\theta_{N}', '\theta_{H}','\theta_{NH_{3}}','\theta_{*}','\theta_{NH}','\theta_{NH_{2}}')
legend boxoff
legend('location','northwest');
toc    