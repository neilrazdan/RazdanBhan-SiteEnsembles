clear all; clc
tic
%% Variables and Parameters
KA = logspace(-1,1,10); % k_ads array
kado = 10^0; % k_des for A2
kro = 10^0; % k_r
kbo = 10^0; % k_ads for B
kbdo = 1; % k_des for B
kadiffo = 10^3; % k_diff for a
kbdiffo = 10^3; % k_diff for b
tspan = [0 10^20]; % time for ode solver
deldel = 10^-4; % differential change for calculation of thermodynamic DoRC
delta = 10^-3; % differential change for calculation of reaction sensitivities
eps = 10^-2; % 1% A and B on surface
%% Arrays for results
S = zeros(length(KA),7)'; r = zeros(1,length(KA)); A = r; O = r; sigma_app = r; TRC_A = r; TRC_B = r; AO = A; AA = A; OO = A; AAA = A; AAO = A; AOA = A; AOO = A; OAO = A; OOO = A;
%% Solvers and Nested Loops
for w = 1:length(KA)
    disp(w)
    for m = 1:1
        k0 = [KA(w) kado kbo kbdo kro kadiffo kbdiffo];
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
                kadiff = K(i,6);
                kbdiff = K(i,7);
                
                %% Sites and Rate Equations
                a = @(x) x(1); b = @(x) x(2);
                aa = @(x) x(3); ab = @(x) x(4); bb = @(x) x(5);
                bbb = @(x) x(6); bba = @(x) x(7);
                abo = @(x) x(8);
                ooo = @(x) x(9); aoo = @(x) x(10);
                aoa = @(x) x(11);
                
                o = @(x) 1 - a(x) - b(x);
                ao = @(x) a(x) - aa(x) - ab(x);
                bo = @(x) b(x) - bb(x) - ab(x);
                oo = @(x) o(x) - ao(x) - bo(x);
                bbo = @(x) bb(x) - bbb(x) - bba(x);
                obo = @(x) bo(x) - bbo(x) - abo(x);
                aba = @(x) ab(x) - bba(x) - abo(x);
                boo = @(x) oo(x) - aoo(x) - ooo(x);
                aob = @(x) ao(x) - aoo(x) - aoa(x);
                bob = @(x) bo(x) - boo(x) - aob(x);
                
                aao = @(x) ao(x);
                aab = @(x) ab(x);
                aaa = @(x) aa(x) - aao(x) - aab(x);
                
                oao = @(x) 0; bao = @(x) 0; bab = @(x) 0;

                % A adsorption
                bboo = @(x) boo(x).*bbo(x)./bo(x);
                oboo = @(x) boo(x).*obo(x)./bo(x);
                aboo = @(x) boo(x).*abo(x)./bo(x);
                oooo = @(x) ooo(x).*ooo(x)./oo(x);
                aooo = @(x) ooo(x).*aoo(x)./oo(x);
                
                ads_a = @(x) 2*ka.*oo(x);
                ads_b = @(x) 0;
                ads_aa = @(x) 2*ka.*(oo(x) + 2*aoo(x));
                ads_ab = @(x) ka*boo(x);
                ads_bb = @(x) 0;
                ads_bbb = @(x) 0;
                ads_bba = @(x) ka*bboo(x);
                ads_abo = @(x) ka.*(oboo(x) - aboo(x));
                ads_ooo = @(x) -2*ka.*(ooo(x) + oooo(x));
                ads_aoo = @(x) ka*(-aoo(x) + oooo(x) - aooo(x));
                ads_aoa = @(x) 2*ka*aooo(x);
                
                % A desorption
                bbaa = @(x) aab(x).*bba(x)./ab(x);
                abaa = @(x) aab(x).*aba(x)./ab(x);
                obaa = @(x) aab(x).*abo(x)./ab(x);
                aaoo = @(x) aoo(x).*aao(x)./ao(x);
                aoaa = @(x) aao(x).*aoa(x)./ao(x);
                
                des_a = @(x) -2*kad.*aa(x);
                des_b = @(x) 0;
                des_aa = @(x) -2*kad.*(aa(x) + 2*aaa(x));
                des_ab = @(x) -kad*aab(x);
                des_bb = @(x) 0;
                des_bbb = @(x) 0;
                des_bba = @(x) -kad*bbaa(x);
                des_abo = @(x) kad.*(abaa(x) - obaa(x));
                des_ooo = @(x) 2*kad.*(aao(x) + aaoo(x));
                des_aoo = @(x) kad*(aaa(x) + aoaa(x) - aaoo(x)); 
                des_aoa = @(x) -2*kad*aoaa(x);
                
                % B adsorption
                b_ads_a = @(x) 0;
                b_ads_b = @(x) kb*o(x);
                b_ads_aa = @(x) 0;
                b_ads_ab = @(x) kb*ao(x);
                b_ads_bb = @(x) 2*kb*bo(x);
                b_ads_bbb = @(x) kb*(bob(x) + 2*bbo(x));
                b_ads_bba = @(x) kb*(aob(x) + abo(x));
                b_ads_abo = @(x) kb.*(aoo(x) - abo(x));
                b_ads_ooo = @(x) -3*kb.*ooo(x);
                b_ads_aoo = @(x) -2*kb.*aoo(x);
                b_ads_aoa = @(x) -kb*aoa(x);
                
                % B desorption
                b_des_a = @(x) 0;
                b_des_b = @(x) -kbd*b(x);
                b_des_aa = @(x) 0;
                b_des_ab = @(x) -kbd*ab(x);
                b_des_bb = @(x) -2*kbd*bb(x);
                b_des_bbb = @(x) -3*kbd*bbb(x);
                b_des_bba = @(x) -2*kbd*bba(x);
                b_des_abo = @(x) kbd.*(bba(x) - abo(x));
                b_des_ooo = @(x) kbd*(obo(x) + 2*boo(x));
                b_des_aoo = @(x) kbd.*(aob(x) + abo(x));
                b_des_aoa = @(x) kbd*aba(x);
                
                % Surface Reaction (aab -> ooo)
                aaab = @(x) aab(x).*aaa(x)./aa(x);
                aaba = @(x) aba(x).*aab(x)./ab(x);
                baab = @(x) aab(x).*aab(x)./aa(x);
                aabb = @(x) aab(x).*bba(x)./ab(x);
                aabbb = @(x) aab(x).*bba(x)./ab(x).*bbb(x)./bb(x);
                aabba = @(x) aab(x).*bba(x)./ab(x).*bba(x)./bb(x);
                bbaab = @(x) bba(x).*aab(x)./ab(x).*aab(x)./aa(x);
                aabo = @(x) abo(x).*aab(x)./ab(x);
                abaab = @(x) aba(x).*aab(x)./ab(x).*aab(x)./aa(x);
                obaab = @(x) abo(x).*aab(x)./ab(x).*aab(x)./aa(x);
                oaab = @(x) aab(x).*aao(x)./aa(x);
                aaboo = @(x) aab(x).*abo(x)./ab(x).*boo(x)./bo(x);
                ooaab = @(x) aoo(x).*aa(x)./ao(x).*aab(x)./aa(x);
                aoaab = @(x) aoa(x).*aao(x)./ao(x).*aab(x)./aa(x);
                aaboa = @(x) aab(x).*abo(x)./ab(x).*aob(x)./bo(x);
                
                rxn_a = @(x) -4*kr*aab(x);
                rxn_b = @(x) -2*kr*aab(x);
                rxn_aa = @(x) -2*kr*(aab(x) + aaab(x));
                rxn_ab = @(x) -kr*(aab(x) + aaba(x) + baab(x));
                rxn_bb = @(x) -2*kr*aabb(x);
                rxn_bbb = @(x) -2*kr*aabbb(x);
                rxn_bba = @(x) -kr*(aabb(x) + aabba(x) + bbaab(x));
                rxn_abo = @(x) kr.*(aabba(x) - aabo(x) + abaab(x) - obaab(x));
                rxn_ooo = @(x) 2*kr.*(aab(x) + aabo(x) + oaab(x) + aaboo(x) + ooaab(x));
                rxn_aoo = @(x) kr.*(aaab(x) + aaba(x) + aoaab(x) + aaboa(x) - ooaab(x));
                rxn_aoa = @(x) -2*kr*aoaab(x);
                
                % Surface Diffusion (aao -> oaa)
                aoaa = @(x) aao(x).*aoa(x)./ao(x);
                aaao = @(x) aao(x).*aaa(x)./aa(x);
                boaa = @(x) aao(x).*aob(x)./ao(x);
                baao = @(x) aao(x).*aab(x)./aa(x);
                aaobb = @(x) aao(x).*aob(x)./ao(x).*bbo(x)./bo(x);
                oaabb = @(x) aao(x).*aab(x)./aa(x).*bba(x)./ab(x);
                oboaa = @(x) aao(x).*aob(x)./ao(x).*obo(x)./bo(x);
                abaao = @(x) aao(x).*aab(x)./aa(x).*aba(x)./ab(x);
                obaao = @(x) aao(x).*aab(x)./aa(x).*abo(x)./ab(x);
                aboaa = @(x) aao(x).*aob(x)./ao(x).*abo(x)./bo(x);
                ooaao = @(x) aao(x).*aao(x)./aa(x).*aoo(x)./ao(x);
                oooaa = @(x) aao(x).*aoo(x)./ao(x).*ooo(x)./oo(x);
                oaao = @(x) aao(x).*aao(x)./aa(x);
                aoaao = @(x) aao(x).*aao(x)./aa(x).*aoa(x)./ao(x);
                aaoo = @(x) aao(x).*aoo(x)./ao(x);
                aooaa = @(x) aao(x).*aoo(x)./ao(x).*aoo(x)./oo(x);
                
                diff_a = @(x) 0;
                diff_b = @(x) 0;
                diff_aa = @(x) 2*kadiff*(aoaa(x) - aaao(x));
                diff_ab = @(x) kadiff*(boaa(x) - baao(x));
                diff_bb = @(x) 0;
                diff_bbb = @(x) 0;
                diff_bba = @(x) kadiff*(aaobb(x) - oaabb(x));
                diff_abo = @(x) kadiff*(oboaa(x) + abaao(x) - obaao(x) - aboaa(x));
                diff_ooo = @(x) 2*kadiff*(ooaao(x) - oooaa(x));
                diff_aoo = @(x) kadiff*(oaao(x) + aoaao(x) + oooaa(x) - aaoo(x) - aooaa(x) - ooaao(x));
                diff_aoa = @(x) 2*kadiff*(aaao(x) + aooaa(x) - aoaa(x) - aoaao(x));
                
                % Surface Diffusion (bo -> ob)
                bbob = @(x) bbo(x).*bob(x)./bo(x);
                bbbo = @(x) bbb(x).*bbo(x)./bb(x);
                abob = @(x) abo(x).*bob(x)./bo(x);
                abbo = @(x) bba(x).*bbo(x)./bb(x);
                oobo = @(x) boo(x).*obo(x)./bo(x);
                ooob = @(x) ooo(x).*boo(x)./oo(x);
                aobo = @(x) aob(x).*obo(x)./bo(x);
                aoob = @(x) aoo(x).*boo(x)./oo(x);
                
                b_diff_a = @(x) 0;
                b_diff_b = @(x) 0;
                b_diff_aa = @(x) 0;
                b_diff_ab = @(x) kbdiff*(aob(x) - abo(x));
                b_diff_bb = @(x) 2*kbdiff*(bob(x) - bbo(x));
                b_diff_bbb = @(x) 2*kbdiff*(bbob(x) - bbbo(x));
                b_diff_bba = @(x) kbdiff*(abob(x) - abbo(x));
                b_diff_abo = @(x) kbdiff*(aob(x) - abo(x) + abbo(x) - abob(x));
                b_diff_ooo = @(x) 2*kbdiff*(oobo(x) - ooob(x));
                b_diff_aoo = @(x) kbdiff*(aobo(x) - aoob(x));
                b_diff_aoa = @(x) 0;
                
                F = @(t,x) [ads_a(x) + des_a(x) + b_ads_a(x) + b_des_a(x) + rxn_a(x) + diff_a(x) + b_diff_a(x);
                    ads_b(x) + des_b(x) + b_ads_b(x) + b_des_b(x) + rxn_b(x) + diff_b(x) + b_diff_b(x);
                    ads_aa(x) + des_aa(x) + b_ads_aa(x) + b_des_aa(x) + rxn_aa(x) + diff_aa(x) + b_diff_aa(x);
                    ads_ab(x) + des_ab(x) + b_ads_ab(x) + b_des_ab(x) + rxn_ab(x) + diff_ab(x) + b_diff_ab(x);
                    ads_bb(x) + des_bb(x) + b_ads_bb(x) + b_des_bb(x) + rxn_bb(x) + diff_bb(x) + b_diff_bb(x);
                    ads_bbb(x) + des_bbb(x) + b_ads_bbb(x) + b_des_bbb(x) + rxn_bbb(x) + diff_bbb(x) + b_diff_bbb(x);
                    ads_bba(x) + des_bba(x) + b_ads_bba(x) + b_des_bba(x) + rxn_bba(x) + diff_bba(x) + b_diff_bba(x);
                    ads_abo(x) + des_abo(x) + b_ads_abo(x) + b_des_abo(x) + rxn_abo(x) + diff_abo(x) + b_diff_abo(x);
                    ads_ooo(x) + des_ooo(x) + b_ads_ooo(x) + b_des_ooo(x) + rxn_ooo(x) + diff_ooo(x) + b_diff_ooo(x);
                    ads_aoo(x) + des_aoo(x) + b_ads_aoo(x) + b_des_aoo(x) + rxn_aoo(x) + diff_aoo(x) + b_diff_aoo(x);
                    ads_aoa(x) + des_aoa(x) + b_ads_aoa(x) + b_des_aoa(x) + rxn_aoa(x) + diff_aoa(x) + b_diff_aoa(x)];
                
                C0 = [eps eps eps^2 eps^2 eps^2 eps^3 eps^3 eps^2*(1-2*eps) (1-2*eps)^3 eps*(1-2*eps)^2 eps^2*(1-2*eps)];                
                
                [T,C] = ode23s(F,tspan,C0);
                t = T;
                a = C(:,1); b = C(:,2);
                aa = C(:,3); ab = C(:,4); bb = C(:,5);
                bbb = C(:,6); bba = C(:,7);
                abo = C(:,8);
                ooo = C(:,9); aoo = C(:,10);
                aoa = C(:,11);
                
                o = 1 - a - b;
                ao = a - aa - ab;
                bo = b - bb - ab;
                oo = o - ao - bo;
                bbo = bb - bbb - bba;
                obo = bo - bbo - abo;
                aba = ab - bba - abo;
                boo = oo - aoo - ooo;
                aob = ao - aoo - aoa;
                bob = bo - boo - aob;
                
                aao = ao;
                aab = ab;
                aaa = aa - aao - aab;
                
                oao = zeros(size(aa)); bao = zeros(size(aa)); bab = zeros(size(aa));
                
                
                aoaa = aao.*aoa./ao;
                aaao = aao.*aaa./aa;
                boaa = aao.*aob./ao;
                baao = aao.*aab./aa;
                aaobb = aao.*aob./ao.*bbo./bo;
                oaabb = aao.*aab./aa.*bba./ab;
                oboaa = aao.*aob./ao.*obo./bo;
                abaao = aao.*aab./aa.*aba./ab;
                obaao = aao.*aab./aa.*abo./ab;
                aboaa = aao.*aob./ao.*abo./bo;
                ooaao = aao.*aao./aa.*aoo./ao;
                oooaa = aao.*aoo./ao.*ooo./oo;
                oaao = aao.*aao./aa;
                aoaao = aao.*aao./aa.*aoa./ao;
                aaoo = aao.*aoo./ao;
                aooaa = aao.*aoo./ao.*aoo./oo;
                
                bbob = bbo.*bob./bo;
                bbbo = bbb.*bbo./bb;
                abob = abo.*bob./bo;
                abbo = bba.*bbo./bb;
                oobo = boo.*obo./bo;
                ooob = ooo.*boo./oo;
                aobo = aob.*obo./bo;
                aoob = aoo.*boo./oo;
                
                z_diff_aa = aaao./aoaa;
                z_diff_ab = baao./boaa;
                z_diff_bba = oaabb./aaobb;
                z_diff_abo = (obaao + aboaa)./(oboaa + abaao);
                z_diff_ooo = oooaa./ooaao;
                z_diff_aoo = (aaoo + aooaa + ooaao)./(oaao + aoaao + oooaa);
                z_diff_aoa = (aoaa + aoaao)./(aaao + aooaa);
                
                z_bdiff_ab = abo./aob;
                z_bdiff_bb = bbo./bob;
                z_bdiff_bbb = bbbo./bbob;
                z_bdiff_bba = abbo./abob;
                z_bdiff_abo = (abo+abob)./(aob+abbo);
                z_bdiff_ooo = ooob./oobo;
                z_bdiff_aoo = aoob./aobo;
                
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
            AAA(m,w) = mean(aaa((end-em):end)); 
            AAO(m,w) = mean(aao((end-em):end)); 
            AOO(m,w) = mean(aoo((end-em):end)); 
            OOO(m,w) = mean(ooo((end-em):end)); 
            AOA(m,w) = mean(aoa((end-em):end)); 
            AAB(m,w) = mean(aab((end-em):end)); 
            
            Zdiff_AA(m,w) = mean(z_diff_aa((end-em):end));
            Zdiff_AB(m,w) = mean(z_diff_ab((end-em):end));
            Zdiff_BBA(m,w) = mean(z_diff_bba((end-em):end));
            Zdiff_ABO(m,w) = mean(z_diff_abo((end-em):end));
            Zdiff_OOO(m,w) = mean(z_diff_ooo((end-em):end));
            Zdiff_AOO(m,w) = mean(z_diff_aoo((end-em):end));
            Zdiff_AOA(m,w) = mean(z_diff_aoa((end-em):end));
            
            Zbdiff_AB(m,w) = mean(z_bdiff_ab((end-em):end));
            Zbdiff_BB(m,w) = mean(z_bdiff_bb((end-em):end));
            Zbdiff_BBB(m,w) = mean(z_bdiff_bbb((end-em):end));
            Zbdiff_BBA(m,w) = mean(z_bdiff_bba((end-em):end));
            Zbdiff_ABO(m,w) = mean(z_bdiff_abo((end-em):end));
            Zbdiff_OOO(m,w) = mean(z_bdiff_ooo((end-em):end));
            Zbdiff_AOO(m,w) = mean(z_bdiff_aoo((end-em):end));
        end
    end
end
%% Extract results
B = 1 - A - O;
AB = A - AO - AA;
BO = O - AO - OO;
BB = B - AB - BO;

% reversibilities
zads = kado.*AA(1,:)./(KA.*OO(1,:));
zads_AA = kado.*(AA(1,:) + 2*AAA(1,:))./(KA.*(OO(1,:) + 2*AOO(1,:)));

% mean-field metrics
mu_AA = AA(1,:)./A(1,:)./A(1,:);
mu_AO = AO(1,:)./A(1,:)./O(1,:);
mu_OO = OO(1,:)./O(1,:)./O(1,:);

mu_AB = AB(1,:)./A(1,:)./B(1,:);
mu_BB = BB(1,:)./B(1,:)./B(1,:);
mu_BO = BO(1,:)./B(1,:)./O(1,:);

mu_AAA = AAA(1,:)./A(1,:)./A(1,:)./A(1,:);
mu_AAB = AAB(1,:)./A(1,:)./A(1,:)./B(1,:);
mu_AAO = AAO(1,:)./A(1,:)./A(1,:)./O(1,:);
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
plot(log10(KA),mu_AB,'r-','linewidth',1.5)
plot(log10(KA),ones(size(mu_AA)),'m--','linewidth',0.5)
xlabel('log_{10}k_{ads}');
ylabel('\mu_{\it aj}');
legend('\mu_{\it aa}','\mu_{\it ao}','\mu_{\it ab}','unity')
legend boxoff
legend('location','northwest');

%% DORC
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),S(1,:)+S(2,:),'k-','linewidth',1.5);
plot(log10(KA),S(3,:)+S(4,:),'r--','linewidth',1.5);
plot(log10(KA),S(5,:),'b.-','linewidth',1.5);
plot(log10(KA),S(6,:),'m:','linewidth',1.5);
plot(log10(KA),S(7,:),'g-','linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('{\it X_{RC,A,ads}}', '{\it X_{RC,B,ads}}', '{\it X_{RC,rxn}}','{\it X_{diff,A}}','{\it X_{diff,B}}');
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

%% Diffusion Reversibilites
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),Zdiff_AA(1,:),'linewidth',1.5);
plot(log10(KA),Zdiff_AB(1,:),'linewidth',1.5);
plot(log10(KA),Zdiff_BBA(1,:),'linewidth',1.5);
plot(log10(KA),Zdiff_ABO(1,:),'linewidth',1.5);
plot(log10(KA),Zdiff_OOO(1,:),'linewidth',1.5);
plot(log10(KA),Zdiff_AOO(1,:),'linewidth',1.5);
plot(log10(KA),Zdiff_AOA(1,:),'linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('aa','ab','bba','abo','ooo','aoo','aoa');
legend boxoff
ylim([0 2]);

figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(log10(KA),Zbdiff_BB(1,:),'linewidth',1.5);
plot(log10(KA),Zbdiff_AB(1,:),'linewidth',1.5);
plot(log10(KA),Zbdiff_BBB(1,:),'linewidth',1.5);
plot(log10(KA),Zbdiff_BBA(1,:),'linewidth',1.5);
plot(log10(KA),Zbdiff_ABO(1,:),'linewidth',1.5);
plot(log10(KA),Zbdiff_OOO(1,:),'linewidth',1.5);
plot(log10(KA),Zbdiff_AOO(1,:),'linewidth',1.5);
xlabel('log_{10}k_{A,ads}');
legend('bb','ab','bbb','bbo','abo','ooo','aoo');
legend boxoff
ylim([0 2]);

toc