clear all; clc
tic
%% Variables and Parameters
KA = logspace(-3,6,50); % k_ads array
kado = 1; % k_des
kro = 10^2/2; % k_r
tspan = [0 10^5]; % time for ode solver
delta = 10^-3; % differential change for calculation of reaction sensitivities
eps = 10^-5; % infinitesmal amount of a on surface at t = 0

% for Section S11.3
EAA = log(1+kro./2./kado); % in units kB*T (positive means repulsion)

% for Section S11.4
% EAA = -2; % in units kB*T (negative means attraction)

% adjust kr accordingly
kro = kro/exp(EAA); 

alpha = 2*(1-exp(-EAA));
%% Arrays for results
S = zeros(length(KA),3)'; r = zeros(2,length(KA)); A = r; O = r; AO = r; AA = r; OO = r;
%% Solvers and Nested Loops
for w = 1:length(KA)
    disp(w)
    for m = 1:1
        k0 = [KA(w) kado kro];
        for j = 1:length(k0)
            K = zeros(2,length(k0));
            K(1,:) = k0; K(2,:) = k0;
            K(2,j) = k0(j)*(1+delta);
            for i = 1:2
                %% Parameters
                ka = K(i,1);
                kad = K(i,2);
                kr = K(i,3);
                
                %% Sites
                a = @(x) x(1);
                o = @(x) 1 - a(x);
                aa = @(x) a(x) - (1 - (1 - 2*alpha.*a(x).*(1-a(x))).^(1/2))./alpha;
                ao = @(x) a(x) - aa(x);
                oo = @(x) o(x) - ao(x);
                phi = @(x) aa(x)./a(x).*exp(EAA) + ao(x)./a(x);

                % A adsorption
                ra_a = @(x) ka.*o(x) - kad*a(x).*phi(x).^4;
                
                % rxn
                rr_a = @(x) -4*kr.*aa(x).*exp(EAA).*phi(x).^6;
                F = @(t,x) [ra_a(x) + rr_a(x)];
                
                C0 = [eps];
                
                [T,C] = ode15s(F,tspan,C0);
                t = T;
                a = C(:,1);
                o = 1 - a;
                aa = a - (1 - (1 - 2*alpha.*a.*(1-a)).^(1/2))./alpha;
                ao = a - aa;
                oo = o - ao;
                
               phi = aa./a.*exp(EAA) + ao./a;

                em = round(length(aa)*0.20);
                R(i,j) = 4*kr.*mean(aa((end-em):end)).*exp(EAA).*mean(phi((end-em):end)).^6;
                R(i,j) = ka.*mean(o((end-em):end)) - kad*mean(a((end-em):end)).*mean(phi((end-em):end)).^4;
                
%                 if mod(w,round(length(KA))/5) == 0 && j == length(k0) && m == 1 && i == 1
%                     figure
%                     set(gca,'fontweight','bold','fontsize',11,'box','on');
%                     hold on
%                     set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
%                     title(num2str(w))
%                     plot(log10(t),log10(aa),'linewidth',1)
%                     plot(log10(t),log10(ao),'linewidth',1)
%                     plot(log10(t),log10(oo),'linewidth',1)
%                     plot(log10(t),log10(a+o),'linewidth',1)
%                     legend('aa','ao','oo','sum')
%                 end
            end
            S(j,w) = (R(2,j) - R(1,j))./(delta*R(1,j));
            r(m,w) = 4*kr.*mean(aa((end-em):end)).*exp(EAA).*mean(phi((end-em):end)).^6;
            A(m,w) = mean(a((end-em):end));
            O(m,w) = mean(o((end-em):end));
            AO(m,w) = mean(ao((end-em):end)); % ao
            AA(m,w) = mean(aa((end-em):end)); % aa
            OO(m,w) = mean(oo((end-em):end)); % oo
        end
    end
end
%% Extract results
PHI = AA(1,:)./A(1,:).*exp(EAA) + AO(1,:)./A(1,:);
kad_eff = kado*PHI.^4;

% reversibilities
zads = kad_eff.*A(1,:)./(KA.*O(1,:));


% mean-field metrics
mu_AA = AA(1,:)./A(1,:)./A(1,:);
mu_AO = AO(1,:)./A(1,:)./O(1,:);
mu_OO = OO(1,:)./O(1,:)./O(1,:);

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
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,ads}}', '{\it X_{RC,r}}');
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
toc

