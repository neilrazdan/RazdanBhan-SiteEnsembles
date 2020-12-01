clear all; clc
tic
%% Variables and Parameters
KA = logspace(-6,3,25); % k_ads array
kado = 10^0; % k_des
kro = 10^2/2; % k_r
kdiffo = 10^4; % k_diff
tspan = [0 10^10]; % time for ode solver
delta = 10^-5; % differential change for calculation of reaction sensitivities
eps = 10^-5; % infinitesmal amount of a on surface at t = 0

% for Section S11.3
EAA = log(1+kro./2./kado); % in units kB*T (positive means repulsion)

% for Section S11.4
EAA = -2; % in units kB*T (negative means attraction)

% adjust kr accordingly
kro = kro/exp(EAA); 
%% Arrays for results
S = zeros(length(KA),4)'; r = zeros(2,length(KA)); A = r; O = r; AO = A; AA = A; OO = A;
%% Solvers and Nested Loops
for w = 1:length(KA)
    disp(w)
    for m = 1:1
        k0 = [KA(w) kado kro kdiffo];
        for j = 1:length(k0)
            K = zeros(2,length(k0));
            K(1,:) = k0; K(2,:) = k0;
            K(2,j) = k0(j)*(1+delta);
            for i = 1:2
                %% Rate Constants
                ka = K(i,1);
                kad = K(i,2);
                kr = K(i,3);
                kdiff = K(i,4);
                
                %% Sites and Rate Equations
                oo = @(x) x(1); aa = @(x) x(2);
                ao = @(x) 1/2 - (aa(x) + oo(x))/2;
                a = @(x) aa(x) + ao(x);
                o = @(x) oo(x) + ao(x);
                
                phi = @(x) aa(x)./a(x).*exp(EAA) + ao(x)./a(x);
        
                % A adsorption
                ads_aa = @(x) 2*ka.*ao(x);
                ads_oo = @(x) -2*ka.*oo(x);
                
                % A desorption
                des_aa = @(x) -2*kad.*aa(x).*exp(EAA).*phi(x).^3;
                des_oo = @(x) 2*kad.*ao(x).*phi(x).^3;
                
                % rxn
                r_aa = @(x) -kr.*aa(x).*(ao(x)./a(x).*exp(EAA).*phi(x).^5 + 7*aa(x)./a(x).*exp(2*EAA).*phi(x).^5);
                r_oo = @(x) kr.*aa(x).*(7*ao(x)./a(x).*exp(EAA).*phi(x).^5 + aa(x)./a(x).*exp(2*EAA).*phi(x).^5);
                
                % diff
                d_aa = @(x) 6*kdiff.*ao(x).*(ao(x)./o(x).*phi(x).^3 - aa(x)./a(x).*exp(EAA).*phi(x).^2);
                d_oo = @(x) 6*kdiff.*ao(x).*(ao(x)./a(x).*phi(x).^2 - oo(x)./o(x).*phi(x).^3);
                
                F = @(t,x) [ads_oo(x) + des_oo(x) + r_oo(x) + d_oo(x); ads_aa(x) + des_aa(x) + r_aa(x) + d_aa(x)];
                
                C0 = [(1-eps)^2 eps^2];
              
                [T,C] = ode15s(F,tspan,C0);
                t = T;
                oo = C(:,1);
                aa = C(:,2);
                ao = 1/2 - (aa+oo)/2;
                a = ao + aa;
                o = ao + oo;
                phi = aa./a.*exp(EAA) + ao./a;

                em = round(length(aa)*0.20);
                R(i,j) = 4*kr.*mean(aa((end-em):end)).*exp(EAA).*mean(phi((end-em):end)).^6;
                
%%%% Histories  
%                 if mod(w,round(length(KA))/5) == 0 && j == length(k0) && m == 1 && i == 1
%                     figure
%                     set(gca,'fontweight','bold','fontsize',11,'box','on');
%                     hold on
%                     set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
%                     title(num2str(w))
%                     plot(log10(t),aa,'linewidth',1)
%                     plot(log10(t),ao,'linewidth',1)
%                     plot(log10(t),oo,'linewidth',1)
%                     plot(log10(t),a+o,'linewidth',1)
%                     legend('aa','ao','oo','sum')
%                 end
            end
            % Results
            S(j,w) = (R(2,j) - R(1,j))./(delta*R(1,j)); % reaction sensitivities
            r(m,w) = 4*kr.*mean(aa((end-em):end)).*exp(EAA).*mean(phi((end-em):end)).^6; % rate
            A(m,w) = mean(a((end-em):end)); % a
            AO(m,w) = mean(ao((end-em):end)); % ao
            AA(m,w) = mean(aa((end-em):end)); % aa
            OO(m,w) = mean(oo((end-em):end)); % oo
            O(m,w) = mean(o((end-em):end)); % o
        end
    end
end
%% Extract results
% Interaction parameters
PHI = AA(1,:)./A(1,:).*exp(EAA) + AO(1,:)./A(1,:);
kad_eff = kado*PHI.^4;

% mean-field metrics
mu_AA = AA(1,:)./A(1,:)./A(1,:);
mu_AO = AO(1,:)./A(1,:)./O(1,:);
mu_OO = OO(1,:)./O(1,:)./O(1,:);

% reversibilities
zads = kad_eff.*A(1,:)./(KA.*O(1,:));
zads_AA = kado.*AA(1,:).*exp(EAA).*PHI(1,:).^3./(KA.*AO(1,:));
zads_OO = kado.*AO(1,:).*PHI(1,:).^3./(KA.*OO(1,:));

% reversibilities
zdiff = AO(1,:)./AO(1,:);
zdiff_AA = (AA(1,:)./A(1,:).*exp(EAA).*PHI.^2)./(AO(1,:)./O(1,:).*PHI.^3);
zdiff_OO = (OO(1,:)./O(1,:).*PHI.^3)./(AO(1,:)./A(1,:).*PHI.^2);

% Quasi-chemical quotient
QC_quot = AA(1,:).*OO(1,:)./(AO(1,:).^2);
QC_K = exp(-EAA);
zQC = QC_quot./QC_K;
%% DoRC sum to unity
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),sum(S,1),'k-','linewidth',1)
ylim([0 2]);
title('DORC check');
%% Mean-field metrics
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),log10(mu_AA),'k--','linewidth',1.5)
plot(log10(KA),log10(mu_AO),'b-.','linewidth',1.5)
plot(log10(KA),log10(mu_OO),'r-','linewidth',1.5)
xlabel('log_{10}k_{ads}');
ylabel('log_{10}\mu_{ij}');
legend('\mu_{aa}','\mu_{ao}','\mu_{oo}')
legend boxoff
legend('location','northwest');
%% Adsorption reversibilities
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),zads,'m:','linewidth',1.5)
plot(log10(KA),zads_AA,'r-','linewidth',1.5)
plot(log10(KA),zads_OO,'b-.','linewidth',1.5)
xlabel('log_{10}k_{ads}');
ylabel('z_{ads,[ij]}');
legend('z_{ads}','z_{ads,[aa]}','z_{ads,[oo]}')
legend boxoff

%% Diffusion reversibilities
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),zdiff,'m:','linewidth',1.5)
plot(log10(KA),zdiff_AA,'r-','linewidth',1.5)
plot(log10(KA),zdiff_OO,'b-.','linewidth',1.5)
xlabel('log_{10}k_{ads}');
ylabel('z_{diff,[ij]}');
legend('z_{diff}','z_{diff,[aa]}','z_{diff,[oo]}')
legend boxoff

%% Reaction Order
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S(1,:),'k-','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('{\it \Phi_{A}}');
legend boxoff

%% Reversibility and sensitivity ratio
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),zads,'k-','linewidth',1.5);
plot(log10(KA),-S(2,:)./S(1,:),'r--','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('{\it z_{ads}}', '{\it -s_{des}/s_{ads}}');
legend boxoff

%% DORC
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S(1,:)+S(2,:),'k-','linewidth',1.5);
plot(log10(KA),S(3,:),'r--','linewidth',1.5);
plot(log10(KA),S(4,:),'b-.','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,ads}}', '{\it X_{RC,r}}', '{\it X_{RC,diff}}');
legend boxoff

%% Quasi-Chemical Approximation
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),QC_quot./QC_K,'k-','linewidth',1.5);
plot(log10(KA),zads,'m:','linewidth',1.5)
xlabel('log_{10}k_{A,ads}');
legend('z_{QC}', 'z_{ads}');
legend boxoff

%% coverage
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),log10(A(1,:)),'k-','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('log_{10} {\it \theta_{a}}');
legend boxoff

toc