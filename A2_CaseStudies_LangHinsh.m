clear all; clc
tic
%% Variables and Parameters
KA = logspace(-6,3,50); % k_ads array
kado = 1; % k_des
kro = 10^2/2; % k_r
tspan = [0 10^10]; % time for ode solver
deldel = 10^-4; % differential change for calculation of thermodynamic DoRC
delta = 10^-4; % differential change for calculation of reaction sensitivities
eps = 10^-10; % infinitesmal amount of a on surface at t = 0
%% Arrays for results
S = zeros(length(KA),3)'; r = zeros(2,length(KA)); A = r; O = r; sigma_app = r; TRC_A = r; r_trc = r;
%% Solvers and Nested Loops
for w = 1:length(KA)
    disp(w)
    for m = 1:2
        k0 = [KA(w) kado*(1+deldel*(m-1)) kro*(1+2*deldel*(m-1))];
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
                
                % A adsorption
                ra_a = @(x) ka.*o(x) - kad*a(x);
                
                % rxn
                rr_a = @(x) -4*kr.*a(x).^2;
                F = @(t,x) [ra_a(x) + rr_a(x)];
                
                C0 = [eps];
                
                [T,C] = ode23s(F,tspan,C0);
                t = T;
                a = C(:,1);
                o = 1 - a;
                
                em = round(length(a)*0.20);
                R(i,j) = 4.*kr.*mean(a((end-em):end)).^2;
                
                if i == 1 && j == 1
                    r_trc(m,w) = 4*kr.*mean(a((end-em):end)).^2;
                end
            end
            S(j,w) = (R(2,j) - R(1,j))./(delta*R(1,j));
            r(m,w) = 4*kr.*mean(a((end-em):end)).^2;
            sigma_app(m,w) = sum(S(1:2,w))+2*S(3,w);
            TRC_A(m,w) = -(S(2,w)+2*S(3,w));
            A(m,w) = mean(a((end-em):end));
            O(m,w) = mean(o((end-em):end));
        end
    end
end
%% Extract results
TRC_A_direct = -(r(2,:) - r(1,:))./(deldel*r(1,:));
zads = kado.*A(1,:)./(KA.*O(1,:));
%% DoRC sum to unity
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),sum(S,1),'k-','linewidth',1)
ylim([0 2]);
title('DORC check');
%% Thermodynamic DoRC
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),-TRC_A_direct,'ko','markersize',6);
plot(log10(KA),-TRC_A(1,:),'rs','markersize',6);
plot(log10(KA),sigma_app(1,:).*A(1,:),'b-.','linewidth',1.5);
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
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,ads}}', '{\it X_{RC,r}}');
legend boxoff
toc

