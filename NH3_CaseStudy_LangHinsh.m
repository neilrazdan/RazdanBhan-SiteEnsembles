clear all; clc
tic
%% Parameters, etc.
Temp = linspace(350+273,550+273,10);
PN2 = 200/4*0.8;
PH2 = 200*3/4*0.8;
PNH3 = 0.2*200;

S = zeros(length(Temp),12)';
norm = 10^6; % just to keep values around unity
% n, nd, h, hd, r1, r1_rev, r2, r2_rev, r3, r3_rev, dd, d
Arow = 1.46*10^13*ones(8,1);
Act = [2.39*10^4; 1.46*10^14; 3.89*10^6; Arow; 2.46*10^3]; % 1/s
Ea = [0.4; 1.2; 0; 1.1; 1.1; 1.1; 1.3; 1.46; 1.2; 0.75; 0.85; 0]; %eV
kB = 8.617333262*10^-5; % eV/K

tspan = [0 10^20];
delta = 10^-3; % differential change for calculation of reaction sensitivities
eps = 10^-2; % all adsorbates start at 1% coverage
%% Output arrays
A = zeros(1,length(Temp)); B = A; C = A; D = A; H = A; O = A;
r = A;
%% Solvers and Nested Loops
for w = 1:length(Temp)
    disp(w)
    T = Temp(w);
    k00 = Act.*exp(-Ea/kB./T)/norm;
    k0 = k00;
    k0(1) = k00(1)*PN2;
    k0(3) = k00(3)*PH2;
    k0(end) = k00(end)*PNH3;
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
                %% Sites
                % a = n, b = nh, c = nh2, d = nh3, h
                d = @(x) x(1); c = @(x) x(2); b = @(x) x(3); a = @(x) x(4); h = @(x) x(5);
                
                o = @(x) 1 - a(x) - b(x) - c(x) - d(x) - h(x); % overall site balance
                oo = @(x) o(x).^2;
                hh = @(x) h(x).^2;
                aa = @(x) a(x).^2;
                ha = @(x) a(x).*h(x);
                hb = @(x) b(x).*h(x);
                hc = @(x) c(x).*h(x);
                bo = @(x) b(x).*o(x);
                co = @(x) c(x).*o(x);
                do = @(x) d(x).*o(x);
                %% N2 adsorption (oo -> aa)
                % d, c, b, a, h
                n_d = @(x) 0; n_c = @(x) 0; n_b = @(x) 0; n_a = @(x) 4*kn*oo(x); n_h = @(x) 0;   
                %% N2 desorption (aa -> oo)
                nd_d = @(x) 0; nd_c = @(x) 0; nd_b = @(x) 0; nd_a = @(x) -4*knd*aa(x); nd_h = @(x) 0;
                %% H2 adsorption (oo -> hh)
                h_d = @(x) 0; h_c = @(x) 0; h_b = @(x) 0; h_a = @(x) 0; h_h = @(x) 4*kh*oo(x);
                %% H2 desorption (hh -> oo)
                hd_d = @(x) 0; hd_c = @(x) 0; hd_b = @(x) 0; hd_a = @(x) 0; hd_h = @(x) -4*khd*hh(x);
                %% N+H (ha -> ob)
                r1_d = @(x) 0; r1_c = @(x) 0; r1_b = @(x) 4*kr1*ha(x); r1_a = @(x) -4*kr1*ha(x); r1_h = @(x) -4*kr1*ha(x);
                %% NH decomp. (ob -> ha)
                r1_rev_d = @(x) 0; r1_rev_c = @(x) 0; r1_rev_b = @(x) -4*kr1_rev*bo(x); r1_rev_a = @(x) 4*kr1_rev*bo(x); r1_rev_h = @(x) 4*kr1_rev*bo(x);
                %% NH+H (hb -> oc)
                r2_d = @(x) 0; r2_c = @(x) 4*kr2*hb(x); r2_b = @(x) -4*kr2*hb(x); r2_a = @(x) 0; r2_h = @(x) -4*kr2*hb(x);
                %% NH2 -> NH+H (oc -> hb)
                r2_rev_d = @(x) 0; r2_rev_c = @(x) -4*kr2_rev*co(x); r2_rev_b = @(x) 4*kr2_rev*co(x); r2_rev_a = @(x) 0; r2_rev_h = @(x) 4*kr2_rev*co(x);
                %% NH2 + H (ch -> do)
                r3_d = @(x) 4*kr3*hc(x); r3_c = @(x) -4*kr3*hc(x); r3_b = @(x) 0; r3_a = @(x) 0; r3_h = @(x) -4*kr3*hc(x);         
                %% NH3 decomp. (do -> ch)
                r3_rev_d = @(x) -4*kr3_rev*do(x); r3_rev_c = @(x) 4*kr3_rev*do(x); r3_rev_b = @(x) 0; r3_rev_a = @(x) 0; r3_rev_h = @(x) 4*kr3_rev*do(x);       
                %% NH3 desorp. (d -> o)
                dd_d = @(x) -kdd.*d(x); dd_c = @(x) 0; dd_b = @(x) 0; dd_a = @(x) 0; dd_h = @(x) 0;
                %% NH3 adsorp. (o -> d)
                d_d = @(x) kd.*o(x); d_c = @(x) 0; d_b = @(x) 0; d_a = @(x) 0; d_h = @(x) 0;
                %% ODE Matrix 
                F = @(t,x) [n_d(x) + nd_d(x) + h_d(x) + hd_d(x) + r1_d(x) + r1_rev_d(x) + r2_d(x) + r2_rev_d(x) + r3_d(x) + r3_rev_d(x) + dd_d(x) + d_d(x);
                    n_c(x) + nd_c(x) + h_c(x) + hd_c(x) + r1_c(x) + r1_rev_c(x) + r2_c(x) + r2_rev_c(x) + r3_c(x) + r3_rev_c(x) + dd_c(x) + d_c(x);
                    n_b(x) + nd_b(x) + h_b(x) + hd_b(x) + r1_b(x) + r1_rev_b(x) + r2_b(x) + r2_rev_b(x) + r3_b(x) + r3_rev_b(x) + dd_b(x) + d_b(x);
                    n_a(x) + nd_a(x) + h_a(x) + hd_a(x) + r1_a(x) + r1_rev_a(x) + r2_a(x) + r2_rev_a(x) + r3_a(x) + r3_rev_a(x) + dd_a(x) + d_a(x);
                    n_h(x) + nd_h(x) + h_h(x) + hd_h(x) + r1_h(x) + r1_rev_h(x) + r2_h(x) + r2_rev_h(x) + r3_h(x) + r3_rev_h(x) + dd_h(x) + d_h(x)];                         
                C00 = ones(1,5)*eps;
                C0 = [C00];
                %% ODE23s
                [T,Y] = ode23s(F,tspan,C0); t = T;
                d = Y(:,1); c = Y(:,2); b = Y(:,3); a = Y(:,4); h = Y(:,5);
               
                o = 1 - a - b - c - d - h; % overall site balance
                oo = o.^2; aa = a.^2;
                hh = h.^2; ha = h.*a; hb = h.*b; hc = h.*c; bo = b.*o; co = c.*o; do = d.*o;
                 
                em = round(length(oo)*0.10); % take last 10% (steady-state)
                
                R(p,j) = kn.*mean(oo(end-em:end)) - knd.*mean(aa(end-em:end)); % rate
                
%                 if mod(w,round(length(P))/10) == 0 && j == length(k0) && m == 1 && p == 1
%                     figure
%                     set(gca,'fontweight','bold','fontsize',11,'box','on');
%                     hold on
%                     set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
%                     title(num2str(w))
%                     plot(log10(t),c,'linewidth',1)
%                     plot(log10(t),h,'linewidth',1)
%                     plot(log10(t),i,'linewidth',1)
%                     plot(log10(t),o,'linewidth',1)
%                     plot(log10(t),f,'linewidth',1)
%                     plot(log10(t),c+o+h+i+f,'linewidth',1)
%                     legend('c','h','i','o','f','sum')
%                 end
            end
            %% Extract coverages, DORC
            S(j,w) = (R(2,j) - R(1,j))./(delta*R(1,j));
            r(m,w) = R(1,1);
            O(m,w) = mean(o(end-em:end));
            A(m,w) = mean(a(end-em:end));
            B(m,w) = mean(b(end-em:end));
            C(m,w) = mean(c(end-em:end));
            D(m,w) = mean(d(end-em:end));
            H(m,w) = mean(h(end-em:end));
        end
    end
end
%% Extract results
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
% N2 adsorption and H2 desorption are the only contributors to E_A,app
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(Temp,S.*Ea*1.60218e-19*6.0221409e+23/1000,'linewidth',1)
xlabel('T / K');
ylabel('E_{A,app} / kJ mol^{-1}');
legend('1','-1','6','-6','2','-2','3','-3','4','-4','5','-5');
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
xlabel('T / K');
legend('X_{RC,N_{2} ads}','X_{RC,H_{2} ads}','X_{RC,r1}','X_{RC,r2}','X_{RC,r3}','X_{RC,NH_{3} des}');
legend boxoff
legend('location','northwest');
ylim([-0.2 1.2])
%% Coverages
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(Temp,A,'linewidth',1.5)
plot(Temp,H,'linewidth',1.5)
plot(Temp,D,'linewidth',1.5)
plot(Temp,O,'linewidth',1.5)
plot(Temp,B,'linewidth',1.5)
plot(Temp,C,'linewidth',1.5)
xlabel('T / K');
legend('\theta_{N}', '\theta_{H}','\theta_{NH_{3}}','\theta_{*}','\theta_{NH}','\theta_{NH_{2}}')
legend boxoff
legend('location','northwest');
toc