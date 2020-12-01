clear all; clc
tic
%% Variables and Parameters
KA = logspace(-1,1,10); % k_ads array
kado = 10^0; % k_des for A2
kro = 10^0; % k_r
kbo = 10^0; % k_ads for B
kbdo = 1; % k_des for B
tspan = [0 10^10]; % time for ode solver
deldel = 10^-4; % differential change for calculation of thermodynamic DoRC
delta = 10^-4; % differential change for calculation of reaction sensitivities
eps = 10^-10; % infinitesmal amount of a and b on surface at t = 0

%% Arrays for results
S = zeros(length(KA),5)'; r = zeros(1,length(KA)); A = r; O = r; sigma_app = r; TRC_A = r; TRC_B = r; AO = A; AA = A; OO = A; AAA = A; AAO = A; AOA = A; AOO = A; OAO = A; OOO = A;
%% Solvers and Nested Loops
for w = 1:length(KA)
    disp(w)
    for m = 1:1
        k0 = [KA(w) kado kbo kbdo kro];
        for j = 1:length(k0)
            K = zeros(2,length(k0));
            K(1,:) = k0; K(2,:) = k0;
            K(2,j) = k0(j)*(1+delta);
            for i = 1:2
                %% Rate Constants
                ka = K(i,1);
                kad = K(i,2);
                kb = K(i,3);
                kbd = K(i,4);
                kr = K(i,5);
                
                %% Sites and Rate Equations
                a = @(x) x(1); b = @(x) x(2);
                
                o = @(x) 1 - a(x) - b(x);
                
                % A adsorption
                
                ads_a = @(x) 2*ka.*o(x).^2;
                ads_b = @(x) 0;
                
                % A desorption
                
                des_a = @(x) -2*kad.*a(x).^2;
                des_b = @(x) 0;
               
                
                % B adsorption
                b_ads_a = @(x) 0;
                b_ads_b = @(x) kb*o(x);
                
                
                % B desorption
                b_des_a = @(x) 0;
                b_des_b = @(x) -kbd*b(x);
               
                % Surface Reaction (aab -> ooo)
                rxn_a = @(x) -4*kr*a(x).^2.*b(x);
                rxn_b = @(x) -2*kr*a(x).^2.*b(x);
               
                
                F = @(t,x) [ads_a(x) + des_a(x) + b_ads_a(x) + b_des_a(x) + rxn_a(x);
                    ads_b(x) + des_b(x) + b_ads_b(x) + b_des_b(x) + rxn_b(x)];
                
                C0 = [eps eps];
                
                [T,C] = ode23s(F,tspan,C0);
                t = T;
                a = C(:,1); b = C(:,2);
                
                o = 1 - a - b;
                aa = a.^2; ab = a.*b; ao = a.*o;
                bb = b.^2; bo = b.*o; oo = o.^2;
                aab = aa.*b;
               
                em = round(length(aab)*0.20);
                R(i,j) = 6*kr.*mean(aab((end-em):end));
                
% %% Histories  
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
            r(m,w) = 6*kr.*mean(aab((end-em):end)); % rate
            A(m,w) = mean(a((end-em):end)); % a
            AO(m,w) = mean(ao((end-em):end)); % ao
            AA(m,w) = mean(aa((end-em):end)); % aa
            OO(m,w) = mean(oo((end-em):end)); % oo
            O(m,w) = mean(o((end-em):end)); % o
        end
    end
end
%% Extract results
B = 1 - A - O;
AB = A - AO - AA;
BO = O - AO - OO;
BB = B - AB - BO;
AAB = AA.*B;

% mean-field metrics
mu_AA = AA(1,:)./A(1,:)./A(1,:);
mu_AO = AO(1,:)./A(1,:)./O(1,:);
mu_OO = OO(1,:)./O(1,:)./O(1,:);
%% DoRC sum to unity
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
hold on
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
plot(log10(KA),sum(S,1),'k-','linewidth',1)
ylim([0 2]);
title('DORC check');

%% DORC
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S(1,:)+S(2,:),'k-','linewidth',1.5);
plot(log10(KA),S(3,:)+S(4,:),'r--','linewidth',1.5);
plot(log10(KA),S(5,:),'b.-','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,A,ads}}', '{\it X_{RC,B,ads}}', '{\it X_{RC,rxn}}');
legend boxoff

%% Reaction Order
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S(1,:),'k-','linewidth',1.5);
plot(log10(KA),S(3,:),'r-','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('{\it \Phi_{A}}','{\it \Phi_{B}}');
legend boxoff

%% Rate
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),log10(AAB(1,:)),'k-','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
ylabel('log_{10} \theta_{\it aab}');
toc