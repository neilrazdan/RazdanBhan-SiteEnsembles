clear all; close all; clc
tic
%% Variables and Parameters
KA = logspace(-6,3,50); % k_ads array
kado = 10^0; % k_des
kro = 10^2/2; % k_r
kdiffo = 10^4; % k_diff
tspan = [0 10^10]; % time for ode solver
deldel = 10^-4; % differential change for calculation of thermodynamic DoRC
delta = 10^-5; % differential change for calculation of reaction sensitivities
eps = 10^-10; % infinitesmal amount of a on surface at t = 0
%% Arrays for results
S = zeros(length(KA),4)'; r = zeros(2,length(KA)); A = r; O = r; sigma_app = r; TRC_A = r; AO = A; AA = A; OO = A;
%% Solvers and Nested Loops
for w = 1:length(KA)
    disp(w)
    for m = 1:2
        k0 = [KA(w) kado*(1+deldel*(m-1)) kro*(1+2*deldel*(m-1)) kdiffo];
        for j = 1:length(k0)
            K = zeros(2,length(k0));
            K(1,:) = k0; K(2,:) = k0;
            K(2,j) = k0(j)*(1+delta);
            for i = 1:2
                %% Rate Constants
                ka = K(i,1);
                kad = K(i,2);
                kr = K(i,3);
                kd = K(i,4);
                %% Sites and Rate Equations
                oo = @(x) x(1); aa = @(x) x(2);
                ao = @(x) 1/2 - (aa(x) + oo(x))/2;
                a = @(x) aa(x) + ao(x);
                o = @(x) oo(x) + ao(x);
        
                % A adsorption
                ads_aa = @(x) 2*ka.*ao(x);
                ads_ao = @(x) ka.*(oo(x)-ao(x));
                ads_oo = @(x) -2*ka.*oo(x);
                
                % A desorption
                des_aa = @(x) -2*kad.*aa(x);
                des_ao = @(x) kad.*(aa(x) - ao(x));
                des_oo = @(x) 2*kad.*ao(x);
                
                % rxn
                r_aa = @(x) -kr.*aa(x).*(1+6.*aa(x)./a(x));
                r_ao = @(x) 3.*kr.*aa(x).*(aa(x)./a(x)-ao(x)./a(x));
                r_oo = @(x) kr.*aa(x).*(1+6.*ao(x)./a(x));
                
                 % diff
                d_aa = @(x) 6*kd.*ao(x).*(ao(x)./o(x) - aa(x)./a(x));
                d_ao = @(x) 3*kd.*ao(x).*(oo(x)./o(x) + aa(x)./a(x) - ao(x)./a(x) - ao(x)./o(x));
                d_oo = @(x) 6*kd.*ao(x).*(ao(x)./a(x) - oo(x)./o(x));
                
                F = @(t,x) [ads_oo(x) + des_oo(x) + r_oo(x) + d_oo(x); ads_aa(x) + des_aa(x) + r_aa(x) + d_aa(x)];
                
                C0 = [(1-eps)^2 eps^2];
              
                [T,C] = ode15s(F,tspan,C0);
                t = T;
                oo = C(:,1);
                aa = C(:,2);
                ao = 1/2 - (aa+oo)/2;
                a = ao + aa;
                o = ao + oo;
                
                em = round(length(aa)*0.20);
                R(i,j) = 4*kr.*mean(aa((end-em):end));
                
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
            sigma_app(m,w) = sum(S(1:2,w))+2*S(3,w) + 2*S(4,w); % apparent sigma
            TRC_A(m,w) = -(S(2,w)+2*S(3,w)); % thermodynamic DoRC of a (surface species)
            r(m,w) = 4*kr.*mean(aa((end-em):end)); % rate
            A(m,w) = mean(a((end-em):end)); % a
            AO(m,w) = mean(ao((end-em):end)); % ao
            AA(m,w) = mean(aa((end-em):end)); % aa
            OO(m,w) = mean(oo((end-em):end)); % oo
            O(m,w) = mean(o((end-em):end)); % o
        end
    end
end
%% Extract results
TRC_A_direct = -(r(2,:) - r(1,:))./(deldel*r(1,:)); % direct calculation of thermodynamic dorc of a

% mean-field metrics
mu_AA = AA(1,:)./A(1,:)./A(1,:);
mu_AO = AO(1,:)./A(1,:)./O(1,:);
mu_OO = OO(1,:)./O(1,:)./O(1,:);

% adsorption reversibilities
zads = kado.*A(1,:)./(KA.*O(1,:));
zads_AA = kado.*AA(1,:)./(KA.*AO(1,:));
zads_OO = kado.*AO(1,:)./(KA.*OO(1,:));
zads_AO = -kado.*(AA(1,:)-AO(1,:))./(KA.*(OO(1,:)-AO(1,:)));

% diffusion reversibilities
zdiff = kdiffo.*AO(1,:)./(kdiffo.*AO(1,:));
zdiff_AA = mu_AA./mu_AO;
zdiff_OO = mu_OO./mu_AO;
zdiff_AO = (mu_AA.*A(1,:)+mu_OO.*O(1,:))./mu_AO;
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
plot(log10(KA),mu_AA,'k--','linewidth',1.5)
plot(log10(KA),mu_AO,'b-.','linewidth',1.5)
plot(log10(KA),mu_OO,'r-','linewidth',1.5)
xlabel('log_{10}k_{ads}');
ylabel('\mu_{ij}');
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
plot(log10(KA),zads_AO,'b-.','linewidth',1.5)
plot(log10(KA),zads_OO,'k--','linewidth',1.5)
xlabel('log_{10}k_{ads}');
ylabel('z_{ads,[ij]}');
legend('z_{ads}','z_{ads,[aa]}','z_{ads,[ao]}','z_{ads,[oo]}')
legend boxoff
%% Diffusion reversibilities
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),zdiff,'m:','linewidth',1.5)
plot(log10(KA),zdiff_AA,'r-','linewidth',1.5)
plot(log10(KA),zdiff_AO,'b-.','linewidth',1.5)
plot(log10(KA),zdiff_OO,'k--','linewidth',1.5)
xlabel('log_{10}k_{ads}');
ylabel('z_{diff,[ij]}');
legend('z_{diff}','z_{diff,[aa]}','z_{diff,[ao]}','z_{diff,[oo]}')
legend boxoff
%% Thermodynamic DoRC
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),-TRC_A_direct,'k-','linewidth',1);
plot(log10(KA),-TRC_A(1,:),'r--','linewidth',1);
plot(log10(KA),sigma_app(1,:).*A(1,:),'b-.','linewidth',1);
xlabel('log_{10}k_{ads}');
legend('{\it -X_{TRC,A}}', '{\it s_{-1} + 2s_{3}}', '{\it \sigma_{*,app}{\theta_{A}}}');
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

%% Reaction Order
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S(1,:),'k-','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('{\it \Phi_{A}}');
legend boxoff

toc