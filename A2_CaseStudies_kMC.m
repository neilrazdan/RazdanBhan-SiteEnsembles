clear all; close all; clc; format long
tic
%% Parameters
KA = 10^(-2)*ones(1,5); % k_ads
d = 256; % lattice size
kAd = 1; % k_des
kr = 10^2/2; % k_r
kdiff = 2*10^-2; % k_diff (multiply by 2 compared to 1st-order)

N = 5*10^6; % number of MC events
t = 1:1:N; % increment MC events
time = zeros(size(t)); % Poisson time

% looping parameters
A = zeros(size(KA)); AA = A; O = A; R = A; ZAA = A; ZAO = A; MUAA = A; AO = A;

% Periodic boundary condition functions
f = @(x) x+1-d*max(0,x+1-d); % right or down periodic boundary condition
g = @(x) x-1-d*min(0,x-2); % left or up periodic boundary condition

% need to calc. z_ads,[aa] and z_ads,[ao]
TON_aa_ads = zeros(1,N);
TON_ao_ads = zeros(1,N);
TON_aa_des = zeros(1,N);
TON_ao_des = zeros(1,N);
TON = zeros(1,N);
%% Big kMC Loop (used to gather repeats)
for w = 1:length(KA)
    %% Setting up
    disp(w)
    toc
    kA = KA(w); % k_ads

    % Lattice matricies
    L = zeros(d,d); % start empty
    Lright = zeros(d,d); Ldown = Lright; % matricies for counting AA reaction pairs
    
    % Initialize arrays for number of a, o, aa, ao
    nA = zeros(1,N); nO = nA; nAA = nA; nAO = nA;
    %% Initialize arrays to keep track of lattice
    [row_A,col_A] = find(L==1); % keep track of A
    [row_O,col_O] = find(L==0); % keep track of O
    
    Lright(:,1:d-1) = L(:,1:d-1)+L(:,2:d); % count horizontal pairs, 0 = oo, 1 = ao, 2 = aa
    Lright(:,d) = L(:,d) + L(:,1); % PBC -- count horizontal pairs, 0 = oo, 1 = ao, 2 = aa
    [row_AA_h,col_AA_h] = find(Lright==2); % keep track of AA pairs
    
    Ldown(1:d-1,:) = L(1:d-1,:) + L(2:d,:); % count vertical pairs, 0 = oo, 1 = ao, 2 = aa
    Ldown(d,:) = L(d,:) + L(1,:); % PBC -- count vertical pairs, 0 = oo, 1 = ao, 2 = aa
    [row_AA_v,col_AA_v] = find(Ldown==2); % keep track of AA pairs
    
    rAA_v = kr.*max(0,length(row_AA_v)); % vertical AA rxn
    rAA_h = kr.*max(0,length(row_AA_h)); % horizontal AA rxn
    rA = kA*max(0,length(row_O)); % A adsorption
    rAd = kAd*max(0,length(row_A)); % A desorption
    
    % diffusion counters and arrays (start empty)
    rAO_down = 0; rAO_up = rAO_down; rAO_right = rAO_down; rAO_left = rAO_down;
    row_AO_up = row_AA_v; row_AO_down = row_AA_v; row_AO_left = row_AA_v; row_AO_right = row_AA_v;
    col_AO_up = row_AA_v; col_AO_down = row_AA_v; col_AO_left = row_AA_v; col_AO_right = row_AA_v;
    
    P = 10^6; % random number range (for picking MC event)
    %% kinetic Monte Carlo
    for i = 1:N
        %     Occasionally display progress
        if any(i/N == [1 2 3 4 5 6 7 8 9 10]/10)
            disp([num2str(i/N*100) '%'])
            disp([num2str(time(i-1))])
            toc
        end
        %% Keeping track of coverages, TON, etc.
        % keep track of TON for ensemble-specific processes
        if i > 1
            TON_aa_ads(i) = TON_aa_ads(i-1);
            TON_ao_ads(i) = TON_ao_ads(i-1);
            TON_aa_des(i) = TON_aa_des(i-1);
            TON_ao_des(i) = TON_ao_des(i-1);
            TON(i) = TON(i-1);
        end
        
        % keep track of coverage, calculate rates
        nA(i) = rAd/kAd; % a = r_des/k_des
        nO(i) = rA/kA; % o = r_ads/k_ads
        nAA(i) = (rAA_h + rAA_v)/kr; % aa = r/k_r
        nAO(i) = (rAO_down + rAO_up + rAO_left + rAO_right)/kdiff;  % ao = r_diff/k_diff
        D = rA + rAd + rAA_v + rAA_h + rAO_down + rAO_up + rAO_left + rAO_right; % total rate
        yAA_v = rAA_v/D; yAA_h = rAA_h/D; yA = rA/D; yAd = rAd/D; % rate process probability
        yAO_down = rAO_down/D; yAO_up = rAO_up/D; yAO_right = rAO_right/D; yAO_left = rAO_left/D; % rate process probability
        dt = -log(rand(1))/D; % Poisson time increment
        time(i) = time(g(i))+dt; % Poisson kMC duration
        
        % Perform MC event
        p = randi([1 P],1); % randomly pick a MC event
        %% Adsorption
        if p <= P*yA
            ads_site = randi([1 round(nO(i))],1);
            site = [row_O(ads_site) col_O(ads_site)]; % coordinates of adsorption site
            
            % If the site above is a-occupied, aa pair is generated
            if isempty(row_A) == 0
                if isempty(find(sum([row_A col_A] == [g(site(1)) site(2)],2)==2,1)) == 0 % "if the site above is a"
                    row_AA_v = [row_AA_v; g(site(1))];
                    col_AA_v = [col_AA_v; site(2)];
                    rAA_v = rAA_v + kr;
                    TON_aa_ads(i) = TON_aa_ads(i) + 1;
                    TON_ao_ads(i) = TON_ao_ads(i) - 1;
                    
                    del_ind = find(sum([row_AO_down col_AO_down] == [g(site(1)) site(2)],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                else
                    row_AO_up = [row_AO_up; site(1)];
                    col_AO_up = [col_AO_up; site(2)];
                    rAO_up = rAO_up + kdiff;
                    
                    TON_ao_ads(i) = TON_ao_ads(i) + 1;
                end
            
            % If the site below is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [f(site(1)) site(2)],2)==2,1)) == 0 % "if the site below is a"
                    row_AA_v = [row_AA_v; site(1)];
                    col_AA_v = [col_AA_v; site(2)];
                    rAA_v = rAA_v + kr;
                    TON_aa_ads(i) = TON_aa_ads(i) + 1;
                    TON_ao_ads(i) = TON_ao_ads(i) - 1;
                    
                    del_ind = find(sum([row_AO_up col_AO_up] == [f(site(1)) site(2)],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                else
                    row_AO_down = [row_AO_down; site(1)];
                    col_AO_down = [col_AO_down; site(2)];
                    rAO_down = rAO_down + kdiff;
                    
                    TON_ao_ads(i) = TON_ao_ads(i) + 1;
                end
            
            % If the site to the right is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [site(1) f(site(2))],2)==2,1)) == 0 % "if the site to the right is a"
                    row_AA_h = [row_AA_h; site(1)];
                    col_AA_h = [col_AA_h; site(2)];
                    rAA_h = rAA_h + kr;
                    TON_aa_ads(i) = TON_aa_ads(i) + 1;
                    TON_ao_ads(i) = TON_ao_ads(i) - 1;
                    
                    del_ind = find(sum([row_AO_left col_AO_left] == [site(1) f(site(2))],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                else
                    row_AO_right = [row_AO_right; site(1)];
                    col_AO_right = [col_AO_right; site(2)];
                    rAO_right = rAO_right + kdiff;
                    
                    TON_ao_ads(i) = TON_ao_ads(i) + 1;
                end
            
            % If the site to the left is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [site(1) g(site(2))],2)==2,1)) == 0 % "if the site to the left is a"
                    row_AA_h = [row_AA_h; site(1)];
                    col_AA_h = [col_AA_h; g(site(2))];
                    rAA_h = rAA_h + kr;
                    TON_aa_ads(i) = TON_aa_ads(i) + 1;
                    TON_ao_ads(i) = TON_ao_ads(i) - 1;
                    
                    del_ind = find(sum([row_AO_right col_AO_right] == [site(1) g(site(2))],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                else
                    row_AO_left = [row_AO_left; site(1)];
                    col_AO_left = [col_AO_left; site(2)];
                    rAO_left = rAO_left + kdiff;
                    TON_ao_ads(i) = TON_ao_ads(i) + 1;
                end
            else
                row_AO_up = [row_AO_up; site(1)];
                col_AO_up = [col_AO_up; site(2)];
                rAO_up = rAO_up + kdiff;
                
                row_AO_down = [row_AO_down; site(1)];
                col_AO_down = [col_AO_down; site(2)];
                rAO_down = rAO_down + kdiff;
                
                row_AO_right = [row_AO_right; site(1)];
                col_AO_right = [col_AO_right; site(2)];
                rAO_right = rAO_right + kdiff;
                
                row_AO_left = [row_AO_left; site(1)];
                col_AO_left = [col_AO_left; site(2)];
                rAO_left = rAO_left + kdiff;
            end
            
            % For all adsorption events
            row_A = [row_A; site(1)]; col_A = [col_A; site(2)]; % update A
            row_O(ads_site) = []; col_O(ads_site) = []; % update O
            rA = rA - kA; % update adsorption rate
            rAd = rAd + kAd; % update desoprtion rate
            
            %% Desorption (Identically, Eley-Rideal reaction)
        elseif p > P*yA && p <= P*(yAd+yA)
            des_site = randi([1 round(nA(i))],1);
            site = [row_A(des_site) col_A(des_site)]; % coordinates of adsorption site
            
            % If the site above is a-occupied, aa pair is annhiliated
            if isempty(row_A) == 0
                if isempty(find(sum([row_A col_A] == [g(site(1)) site(2)],2)==2,1)) == 0 % "if the site above is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [g(site(1)) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    TON_aa_des(i) = TON_aa_des(i) + 1;
                    TON_ao_des(i) = TON_ao_des(i) - 1;
                    
                    row_AO_down = [row_AO_down; g(site(1))];
                    col_AO_down = [col_AO_down; site(2)];
                    rAO_down = rAO_down + kdiff;
                else
                    del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                    
                    TON_ao_des(i) = TON_ao_des(i) + 1;
                end
            
            % If the site below is a-occupied, aa pair is annhiliated
                if isempty(find(sum([row_A col_A] == [f(site(1)) site(2)],2)==2,1)) == 0 % "if the site below is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [site(1) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    TON_aa_des(i) = TON_aa_des(i) + 1;
                    TON_ao_des(i) = TON_ao_des(i) - 1;
                    
                    row_AO_up = [row_AO_up; f(site(1))];
                    col_AO_up = [col_AO_up; site(2)];
                    rAO_up = rAO_up + kdiff;
                else
                    del_ind = find(sum([row_AO_down col_AO_down] == [site(1) site(2)],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                    
                    TON_ao_des(i) = TON_ao_des(i) + 1;
                end
            
            % If the site to the right is a-occupied, aa pair is annhiliated
                if isempty(find(sum([row_A col_A] == [site(1) f(site(2))],2)==2,1)) == 0 % "if the site to the right is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) site(2)],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    TON_aa_des(i) = TON_aa_des(i) + 1;
                    TON_ao_des(i) = TON_ao_des(i) - 1;
                    
                    row_AO_left = [row_AO_left; site(1)];
                    col_AO_left = [col_AO_left; f(site(2))];
                    rAO_left = rAO_left + kdiff;
                else
                    del_ind = find(sum([row_AO_right col_AO_right] == [site(1) site(2)],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                    
                    TON_ao_des(i) = TON_ao_des(i) + 1;
                end
            
            % If the site to the left is a-occupied, aa pair is annhiliated
                if isempty(find(sum([row_A col_A] == [site(1) g(site(2))],2)==2,1)) == 0 % "if the site to the left is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) g(site(2))],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    TON_aa_des(i) = TON_aa_des(i) + 1;
                    TON_ao_des(i) = TON_ao_des(i) - 1;
                    
                    row_AO_right = [row_AO_right; site(1)];
                    col_AO_right = [col_AO_right; g(site(2))];
                    rAO_right = rAO_right + kdiff;
                else
                    del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                    
                    TON_ao_des(i) = TON_ao_des(i) + 1;
                end
            else
                del_ind = find(sum([row_AO_right col_AO_right] == [site(1) site(2)],2)==2,1);
                row_AO_right(del_ind) = [];
                col_AO_right(del_ind) = [];
                rAO_right = rAO_right - kdiff;
                
                del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                row_AO_left(del_ind) = [];
                col_AO_left(del_ind) = [];
                rAO_left = rAO_left - kdiff;
                
                del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                row_AO_up(del_ind) = [];
                col_AO_up(del_ind) = [];
                rAO_up = rAO_up - kdiff;
                
                del_ind = find(sum([row_AO_down col_AO_down] == [site(1) site(2)],2)==2,1);
                row_AO_down(del_ind) = [];
                col_AO_down(del_ind) = [];
                rAO_down = rAO_down - kdiff;
                
                TON_ao_des(i) = TON_ao_des(i) + 4;
            end
            
            % For every desorption event
            row_O = [row_O; site(1)]; col_O = [col_O; site(2)]; % update O
            row_A(des_site) = []; col_A(des_site) = []; % update A
            rA = rA + kA; % update adsorption rate
            rAd = rAd - kAd; % update desoprtion rate
            %% Horizontal Reaction
        elseif p > P*(yAd+yA) && p <= P*(yAd+yA+yAA_h)
            r_h_site = randi([1 round(rAA_h/kr)],1); % random horizontal site
            site = [row_AA_h(r_h_site) col_AA_h(r_h_site)]; % coordinates of rxn site pair
            
            % If the site above the LHS-a in aa pair is a-occupied, aa pair
            % is annihilated
            if isempty(row_A) == 0
                if isempty(find(sum([row_A col_A] == [g(site(1)) site(2)],2)==2,1)) == 0 % "if the site above is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [g(site(1)) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    row_AO_down = [row_AO_down; g(site(1))];
                    col_AO_down = [col_AO_down; site(2)];
                    rAO_down = rAO_down + kdiff;
                else
                    del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                end
            
            % If the site below the LHS-a in aa pair is a-occupied, aa pair
            % is annihilated
                if isempty(find(sum([row_A col_A] == [f(site(1)) site(2)],2)==2,1)) == 0 % "if the site below is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [site(1) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    row_AO_up = [row_AO_up; f(site(1))];
                    col_AO_up = [col_AO_up; site(2)];
                    rAO_up = rAO_up + kdiff;
                else
                    del_ind = find(sum([row_AO_down col_AO_down] == [site(1) site(2)],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                end
            
            % If the site above the RHS-a in aa pair is a-occupied, aa pair
            % is annihilated
                if isempty(find(sum([row_A col_A] == [g(site(1)) f(site(2))],2)==2,1)) == 0 % "if the site R and Up is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [g(site(1)) f(site(2))],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    row_AO_down = [row_AO_down; g(site(1))];
                    col_AO_down = [col_AO_down; f(site(2))];
                    rAO_down = rAO_down + kdiff;
                else
                    del_ind = find(sum([row_AO_up col_AO_up] == [site(1) f(site(2))],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                end
            
            % If the site below the RHS-a in aa pair is a-occupied, aa pair
            % is annihilated
                if isempty(find(sum([row_A col_A] == [f(site(1)) f(site(2))],2)==2,1)) == 0 % "if the site R and Down is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [site(1) f(site(2))],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    row_AO_up = [row_AO_up; f(site(1))];
                    col_AO_up = [col_AO_up; f(site(2))];
                    rAO_up = rAO_up + kdiff;
                else
                    del_ind = find(sum([row_AO_down col_AO_down] == [site(1) f(site(2))],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                end
            
            % If the site to the left of the LHS-a in aa pair is a-occupied, aa pair
            % is annihilated
                if isempty(find(sum([row_A col_A] == [site(1) g(site(2))],2)==2,1)) == 0 % "if the site to the left is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) g(site(2))],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    row_AO_right = [row_AO_right; site(1)];
                    col_AO_right = [col_AO_right; g(site(2))];
                    rAO_right = rAO_right + kdiff;
                else
                    del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                end
            
            % If the site to the right of the RHS-a in aa pair is a-occupied, aa pair
            % is annihilated
                if isempty(find(sum([row_A col_A] == [site(1) f(f(site(2)))],2)==2,1)) == 0 % "if the site 2 to the right is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) f(site(2))],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    row_AO_left = [row_AO_left; site(1)];
                    col_AO_left = [col_AO_left; f(f(site(2)))];
                    rAO_left = rAO_left + kdiff;
                else
                    del_ind = find(sum([row_AO_right col_AO_right] == [site(1) f(site(2))],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                end
            else
                del_ind = find(sum([row_AO_down col_AO_down] == [site(1) f(site(2))],2)==2,1);
                row_AO_down(del_ind) = [];
                col_AO_down(del_ind) = [];
                rAO_down = rAO_down - kdiff;
                
                del_ind = find(sum([row_AO_up col_AO_up] == [site(1) f(site(2))],2)==2,1);
                row_AO_up(del_ind) = [];
                col_AO_up(del_ind) = [];
                rAO_up = rAO_up - kdiff;
                
                del_ind = find(sum([row_AO_down col_AO_down] == [site(1) site(2)],2)==2,1);
                row_AO_down(del_ind) = [];
                col_AO_down(del_ind) = [];
                rAO_down = rAO_down - kdiff;
                
                del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                row_AO_up(del_ind) = [];
                col_AO_up(del_ind) = [];
                rAO_up = rAO_up - kdiff;
                
                del_ind = find(sum([row_AO_right col_AO_right] == [site(1) f(site(2))],2)==2,1);
                row_AO_right(del_ind) = [];
                col_AO_right(del_ind) = [];
                rAO_right = rAO_right - kdiff;
                
                del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                row_AO_left(del_ind) = [];
                col_AO_left(del_ind) = [];
                rAO_left = rAO_left - kdiff;
            end
            
            % add two o
            row_O = [row_O; site(1); site(1)]; col_O = [col_O; site(2); f(site(2))]; % update 2 O's
            
            % delete both a
            del_ind = find(sum([row_A col_A] == [site(1) site(2)],2)==2,1);
            row_A(del_ind) = []; col_A(del_ind) = []; % update 1st A in pair
            del_ind = find(sum([row_A col_A] == [site(1) f(site(2))],2)==2,1);
            row_A(del_ind) = []; col_A(del_ind) = []; % update 2nd A in pair
            
            % delete aa pair
            del_ind = find(sum([row_AA_h col_AA_h] == [site(1) site(2)],2)==2,1);
            row_AA_h(del_ind) = []; col_AA_h(del_ind) = [];
            
            rA = rA + 2*kA; % update adsorption rate
            rAd = rAd - 2*kAd; % update desoprtion rate
            rAA_h = rAA_h - kr; % update horizontal rxn rate
            TON(i) = TON(i) + 1;
            %% Vertical Reaction
        elseif p > P*(yAd+yA+yAA_h) && p <= P*(yAd+yA+yAA_h+yAA_v)
            r_v_site = randi([1 round(rAA_v/kr)],1); % random vertical site
            site = [row_AA_v(r_v_site) col_AA_v(r_v_site)]; % coordinates of rxn site pair
            
            % "if the site above is also a"
            if isempty(row_A) == 0
                if isempty(find(sum([row_A col_A] == [g(site(1)) site(2)],2)==2,1)) == 0 % "if the site above is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [g(site(1)) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    row_AO_down = [row_AO_down; g(site(1))];
                    col_AO_down = [col_AO_down; site(2)];
                    rAO_down = rAO_down + kdiff;
                else
                    del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                end
            
            % "if the site two below is also a"
                if isempty(find(sum([row_A col_A] == [f(f(site(1))) site(2)],2)==2,1)) == 0 % "if the site two below also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [f(site(1)) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    row_AO_up = [row_AO_up; f(f(site(1)))];
                    col_AO_up = [col_AO_up; site(2)];
                    rAO_up = rAO_up + kdiff;
                else
                    del_ind = find(sum([row_AO_down col_AO_down] == [f(site(1)) site(2)],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                end
            
            % "if the site to the right is also a"
                if isempty(find(sum([row_A col_A] == [site(1) f(site(2))],2)==2,1)) == 0 % "if the site to the right is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) site(2)],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    row_AO_left = [row_AO_left; site(1)];
                    col_AO_left = [col_AO_left; f(site(2))];
                    rAO_left = rAO_left + kdiff;
                else
                    del_ind = find(sum([row_AO_right col_AO_right] == [site(1) site(2)],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                end
            
            % "if the site to the left is also a"
                if isempty(find(sum([row_A col_A] == [site(1) g(site(2))],2)==2,1)) == 0 % "if the site to the left is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) g(site(2))],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    row_AO_right = [row_AO_right; site(1)];
                    col_AO_right = [col_AO_right; g(site(2))];
                    rAO_right = rAO_right + kdiff;
                else
                    del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                end
            
            % "if the site down and to the right is also a"
                if isempty(find(sum([row_A col_A] == [f(site(1)) f(site(2))],2)==2,1)) == 0 % "if the site down and to the right is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [f(site(1)) site(2)],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    row_AO_left = [row_AO_left; f(site(1))];
                    col_AO_left = [col_AO_left; f(site(2))];
                    rAO_left = rAO_left + kdiff;
                else
                    del_ind = find(sum([row_AO_right col_AO_right] == [f(site(1)) site(2)],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                end
            
            % "if the site down and to the left is also a"
                if isempty(find(sum([row_A col_A] == [f(site(1)) g(site(2))],2)==2,1)) == 0 % "if the site down and to the left is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [f(site(1)) g(site(2))],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    
                    row_AO_right = [row_AO_right; f(site(1))];
                    col_AO_right = [col_AO_right; g(site(2))];
                    rAO_right = rAO_right + kdiff;
                else
                    del_ind = find(sum([row_AO_left col_AO_left] == [f(site(1)) site(2)],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                end
            else
                del_ind = find(sum([row_AO_left col_AO_left] == [f(site(1)) site(2)],2)==2,1);
                row_AO_left(del_ind) = [];
                col_AO_left(del_ind) = [];
                rAO_left = rAO_left - kdiff;
                
                del_ind = find(sum([row_AO_right col_AO_right] == [f(site(1)) site(2)],2)==2,1);
                row_AO_right(del_ind) = [];
                col_AO_right(del_ind) = [];
                rAO_right = rAO_right - kdiff;
                
                del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                row_AO_left(del_ind) = [];
                col_AO_left(del_ind) = [];
                rAO_left = rAO_left - kdiff;
                
                del_ind = find(sum([row_AO_right col_AO_right] == [site(1) site(2)],2)==2,1);
                row_AO_right(del_ind) = [];
                col_AO_right(del_ind) = [];
                rAO_right = rAO_right - kdiff;
                
                del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                row_AO_up(del_ind) = [];
                col_AO_up(del_ind) = [];
                rAO_up = rAO_up - kdiff;
                
                del_ind = find(sum([row_AO_down col_AO_down] == [f(site(1)) site(2)],2)==2,1);
                row_AO_down(del_ind) = [];
                col_AO_down(del_ind) = [];
                rAO_down = rAO_down - kdiff;
            end
            
            % add two o
            row_O = [row_O; site(1); f(site(1))]; col_O = [col_O; site(2); site(2)]; % update 2 O's
            
            % delete both a
            del_ind = find(sum([row_A col_A] == [site(1) site(2)],2)==2,1);
            row_A(del_ind) = []; col_A(del_ind) = []; % update 1st A in pair
            del_ind = find(sum([row_A col_A] == [f(site(1)) site(2)],2)==2,1);
            row_A(del_ind) = []; col_A(del_ind) = []; % update 2nd A in pair
            
            % delete aa pair
            del_ind = find(sum([row_AA_v col_AA_v] == [site(1) site(2)],2)==2,1);
            row_AA_v(del_ind) = []; col_AA_v(del_ind) = [];
            
            rA = rA + 2*kA; % update adsorption rate
            rAd = rAd - 2*kAd; % update desoprtion rate
            rAA_v = rAA_v - kr; % update vertical rxn rate
            TON(i) = TON(i) + 1;
            %% Right Diffusion
        elseif p > P*(yAd+yA+yAA_h+yAA_v) && p <= P*(yAd+yA+yAA_h+yAA_v+yAO_right)
            right_site = randi([1 round(rAO_right/kdiff)],1); % random site
            site = [row_AO_right(right_site) col_AO_right(right_site)]; % coordinates of site pair
             
            if isempty(row_A) == 0
                % "if the site above is also a"
                if isempty(find(sum([row_A col_A] == [g(site(1)) site(2)],2)==2,1)) == 0 % "if the site above is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [g(site(1)) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    
                    row_AO_down = [row_AO_down; g(site(1))];
                    col_AO_down = [col_AO_down; site(2)];
                    rAO_down = rAO_down + kdiff;
                else
                    del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                end
                
                % If the site below is a-occupied, aa pair is annhiliated
                if isempty(find(sum([row_A col_A] == [f(site(1)) site(2)],2)==2,1)) == 0 % "if the site below is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [site(1) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    
                    row_AO_up = [row_AO_up; f(site(1))];
                    col_AO_up = [col_AO_up; site(2)];
                    rAO_up = rAO_up + kdiff;
                else
                    del_ind = find(sum([row_AO_down col_AO_down] == [site(1) site(2)],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                    
                end
                
                % If the site to the left is a-occupied, aa pair is annhiliated
                if isempty(find(sum([row_A col_A] == [site(1) g(site(2))],2)==2,1)) == 0 % "if the site to the left is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) g(site(2))],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    
                    row_AO_right = [row_AO_right; site(1)];
                    col_AO_right = [col_AO_right; g(site(2))];
                    rAO_right = rAO_right + kdiff;
                    
                else
                    del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                end
                
                % If the site above the diffusion acceptor is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [g(site(1)) f(site(2))],2)==2,1)) == 0 % "if the site above is a"
                    row_AA_v = [row_AA_v; g(site(1))];
                    col_AA_v = [col_AA_v; f(site(2))];
                    rAA_v = rAA_v + kr;
                    
                    
                    del_ind = find(sum([row_AO_down col_AO_down] == [g(site(1)) f(site(2))],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                else
                    row_AO_up = [row_AO_up; site(1)];
                    col_AO_up = [col_AO_up; f(site(2))];
                    rAO_up = rAO_up + kdiff;
                    
                end
                
                % If the site below the diff accept. is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [f(site(1)) f(site(2))],2)==2,1)) == 0 % "if the site below is a"
                    row_AA_v = [row_AA_v; site(1)];
                    col_AA_v = [col_AA_v; f(site(2))];
                    rAA_v = rAA_v + kr;
                    
                    
                    del_ind = find(sum([row_AO_up col_AO_up] == [f(site(1)) f(site(2))],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                else
                    row_AO_down = [row_AO_down; site(1)];
                    col_AO_down = [col_AO_down; f(site(2))];
                    rAO_down = rAO_down + kdiff;
                end
                
                % If the site to the right of the diff. accept. is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [site(1) f(f(site(2)))],2)==2,1)) == 0 % "if the site to the right is a"
                    row_AA_h = [row_AA_h; site(1)];
                    col_AA_h = [col_AA_h; f(site(2))];
                    rAA_h = rAA_h + kr;
                    
                    
                    del_ind = find(sum([row_AO_left col_AO_left] == [site(1) f(f(site(2)))],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                else
                    row_AO_right = [row_AO_right; site(1)];
                    col_AO_right = [col_AO_right; f(site(2))];
                    rAO_right = rAO_right + kdiff;
                end
            end
            
            % delete a and o
            del_ind = find(sum([row_A col_A] == [site(1) site(2)],2)==2,1);
            row_A(del_ind) = []; col_A(del_ind) = []; 
            del_ind = find(sum([row_O col_O] == [site(1) f(site(2))],2)==2,1);
            row_O(del_ind) = []; col_O(del_ind) = []; 
            
            % add a and o
            row_O = [row_O; site(1)]; col_O = [col_O; site(2)]; % update O
            row_A = [row_A; site(1)]; col_A = [col_A; f(site(2))]; % update A
            
            % delete ao pair
            del_ind = find(sum([row_AO_right col_AO_right] == [site(1) site(2)],2)==2,1);
            row_AO_right(del_ind) = []; col_AO_right(del_ind) = [];
            rAO_right = rAO_right - kdiff;
            
            % make ao pair
            row_AO_left = [row_AO_left; site(1)]; col_AO_left = [col_AO_left; f(site(2))]; % update AO
            rAO_left = rAO_left + kdiff;
        %% Left Diffusion
        elseif p > P*(yAd+yA+yAA_h+yAA_v+yAO_right) && p <= P*(yAd+yA+yAA_h+yAA_v+yAO_right+yAO_left)
            left_site = randi([1 round(rAO_left/kdiff)],1); % random site
            site = [row_AO_left(left_site) col_AO_left(left_site)]; % coordinates of site pair
             
            if isempty(row_A) == 0
                % "if the site above is also a"
                if isempty(find(sum([row_A col_A] == [g(site(1)) site(2)],2)==2,1)) == 0 % "if the site above is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [g(site(1)) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    
                    row_AO_down = [row_AO_down; g(site(1))];
                    col_AO_down = [col_AO_down; site(2)];
                    rAO_down = rAO_down + kdiff;
                else
                    del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                end
                
                % "if the site below is also a"
               if isempty(find(sum([row_A col_A] == [f(site(1)) site(2)],2)==2,1)) == 0 % "if the site below is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [site(1) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    
                    row_AO_up = [row_AO_up; f(site(1))];
                    col_AO_up = [col_AO_up; site(2)];
                    rAO_up = rAO_up + kdiff;
                else
                    del_ind = find(sum([row_AO_down col_AO_down] == [site(1) site(2)],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                end
                
                % "if the site to the right is also a"
                if isempty(find(sum([row_A col_A] == [site(1) f(site(2))],2)==2,1)) == 0 % "if the site to the right is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) site(2)],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    
                    row_AO_left = [row_AO_left; site(1)];
                    col_AO_left = [col_AO_left; f(site(2))];
                    rAO_left = rAO_left + kdiff;
                else
                    del_ind = find(sum([row_AO_right col_AO_right] == [site(1) site(2)],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                end
                
                % if the site above the diff. acceptor is also a
                if isempty(find(sum([row_A col_A] == [g(site(1)) g(site(2))],2)==2,1)) == 0 
                    row_AA_v = [row_AA_v; g(site(1))];
                    col_AA_v = [col_AA_v; g(site(2))];
                    rAA_v = rAA_v + kr;
                    
                    
                    del_ind = find(sum([row_AO_down col_AO_down] == [g(site(1)) g(site(2))],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                else
                    row_AO_up = [row_AO_up; site(1)];
                    col_AO_up = [col_AO_up; g(site(2))];
                    rAO_up = rAO_up + kdiff;
                end
                
                % if the site below the diff. acceptor is also a
                if isempty(find(sum([row_A col_A] == [f(site(1)) g(site(2))],2)==2,1)) == 0 % "if the site below is a"
                    row_AA_v = [row_AA_v; site(1)];
                    col_AA_v = [col_AA_v; g(site(2))];
                    rAA_v = rAA_v + kr;
                    
                     
                    del_ind = find(sum([row_AO_up col_AO_up] == [f(site(1)) g(site(2))],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                else
                    row_AO_down = [row_AO_down; site(1)];
                    col_AO_down = [col_AO_down; g(site(2))];
                    rAO_down = rAO_down + kdiff;
                end
                
                % if the site to the left of the diff. acceptor is also a
                if isempty(find(sum([row_A col_A] == [site(1) g(g(site(2)))],2)==2,1)) == 0 % "if the site to the left is a"
                    row_AA_h = [row_AA_h; site(1)];
                    col_AA_h = [col_AA_h; g(g(site(2)))];
                    rAA_h = rAA_h + kr;
                    
                    
                    del_ind = find(sum([row_AO_right col_AO_right] == [site(1) g(g(site(2)))],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                else
                    row_AO_left = [row_AO_left; site(1)];
                    col_AO_left = [col_AO_left; g(site(2))];
                    rAO_left = rAO_left + kdiff;
                end
            end
            
            % delete a and o
            del_ind = find(sum([row_A col_A] == [site(1) site(2)],2)==2,1);
            row_A(del_ind) = []; col_A(del_ind) = [];
            del_ind = find(sum([row_O col_O] == [site(1) g(site(2))],2)==2,1);
            row_O(del_ind) = []; col_O(del_ind) = [];
            
            % add a and o
            row_O = [row_O; site(1)]; col_O = [col_O; site(2)]; % update O
            row_A = [row_A; site(1)]; col_A = [col_A; g(site(2))]; % update A
            
            % delete ao pair
            del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
            row_AO_left(del_ind) = []; col_AO_left(del_ind) = [];
            rAO_left = rAO_left - kdiff;
            
            % make ao pair
            row_AO_right = [row_AO_right; site(1)]; col_AO_right = [col_AO_right; g(site(2))]; % update AO
            rAO_right = rAO_right + kdiff;
            %% Up Diffusion
        elseif p > P*(yAd+yA+yAA_h+yAA_v+yAO_right+yAO_left) && p <= P*(yAd+yA+yAA_h+yAA_v+yAO_right+yAO_left+yAO_up)
            up_site = randi([1 round(rAO_up/kdiff)],1); % random site
            site = [row_AO_up(up_site) col_AO_up(up_site)]; % coordinates of site pair
            
            if isempty(row_A) == 0
                % "if the site below is also a"
                if isempty(find(sum([row_A col_A] == [f(site(1)) site(2)],2)==2,1)) == 0 % "if the site below is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [site(1) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                    
                    
                    row_AO_up = [row_AO_up; f(site(1))];
                    col_AO_up = [col_AO_up; site(2)];
                    rAO_up = rAO_up + kdiff;
                else
                    del_ind = find(sum([row_AO_down col_AO_down] == [site(1) site(2)],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                end
                
                % "if the site to the right is also a"
                if isempty(find(sum([row_A col_A] == [site(1) f(site(2))],2)==2,1)) == 0 % "if the site to the right is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) site(2)],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    
                    row_AO_left = [row_AO_left; site(1)];
                    col_AO_left = [col_AO_left; f(site(2))];
                    rAO_left = rAO_left + kdiff;
                else
                    del_ind = find(sum([row_AO_right col_AO_right] == [site(1) site(2)],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                end
                
                % If the site to the left is a-occupied, aa pair is annhiliated
                if isempty(find(sum([row_A col_A] == [site(1) g(site(2))],2)==2,1)) == 0 % "if the site to the left is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) g(site(2))],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    
                    row_AO_right = [row_AO_right; site(1)];
                    col_AO_right = [col_AO_right; g(site(2))];
                    rAO_right = rAO_right + kdiff;
                    
                else
                    del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                end
                
                % If the site above the diffusion acceptor is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [g(g(site(1))) site(2)],2)==2,1)) == 0 % "if the site above is a"
                    row_AA_v = [row_AA_v; g(g(site(1)))];
                    col_AA_v = [col_AA_v; site(2)];
                    rAA_v = rAA_v + kr;
                    
                    
                    del_ind = find(sum([row_AO_down col_AO_down] == [g(g(site(1))) site(2)],2)==2,1);
                    row_AO_down(del_ind) = [];
                    col_AO_down(del_ind) = [];
                    rAO_down = rAO_down - kdiff;
                else
                    row_AO_up = [row_AO_up; g(site(1))];
                    col_AO_up = [col_AO_up; site(2)];
                    rAO_up = rAO_up + kdiff;
                end
                
                % If the site to the right of the diff. accept. is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [g(site(1)) f(site(2))],2)==2,1)) == 0 % "if the site to the right is a"
                    row_AA_h = [row_AA_h; g(site(1))];
                    col_AA_h = [col_AA_h; site(2)];
                    rAA_h = rAA_h + kr;
                     
                     
                    del_ind = find(sum([row_AO_left col_AO_left] == [g(site(1)) f(site(2))],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                else
                    row_AO_right = [row_AO_right; g(site(1))];
                    col_AO_right = [col_AO_right; site(2)];
                    rAO_right = rAO_right + kdiff;
                end
                
                % If the site to the left of the diff. accept. is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [g(site(1)) g(site(2))],2)==2,1)) == 0 % "if the site to the right is a"
                    row_AA_h = [row_AA_h; g(site(1))];
                    col_AA_h = [col_AA_h; g(site(2))];
                    rAA_h = rAA_h + kr;
                     
                     
                    del_ind = find(sum([row_AO_right col_AO_right] == [g(site(1)) g(site(2))],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                else
                    row_AO_left = [row_AO_left; g(site(1))];
                    col_AO_left = [col_AO_left; site(2)];
                    rAO_left = rAO_left + kdiff;
                end
            end
            
            % delete a and o
            del_ind = find(sum([row_A col_A] == [site(1) site(2)],2)==2,1);
            row_A(del_ind) = []; col_A(del_ind) = [];
            del_ind = find(sum([row_O col_O] == [g(site(1)) site(2)],2)==2,1);
            row_O(del_ind) = []; col_O(del_ind) = [];
            
            % add a and o
            row_O = [row_O; site(1)]; col_O = [col_O; site(2)]; % update O
            row_A = [row_A; g(site(1))]; col_A = [col_A; site(2)]; % update A
            
            % delete ao pair
            del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
            row_AO_up(del_ind) = []; col_AO_up(del_ind) = [];
            rAO_up = rAO_up - kdiff;
            
            % make ao pair
            row_AO_down = [row_AO_down; g(site(1))]; col_AO_down = [col_AO_down; site(2)]; % update AO
            rAO_down = rAO_down + kdiff;
            %% Down Diffusion
        elseif p > P*(yAd+yA+yAA_h+yAA_v+yAO_right+yAO_left+yAO_up) && p <= P*(yAd+yA+yAA_h+yAA_v+yAO_right+yAO_left+yAO_up+yAO_down)
            down_site = randi([1 round(rAO_down/kdiff)],1); % random site
            site = [row_AO_down(down_site) col_AO_down(down_site)]; % coordinates of site pair
             
            if isempty(row_A) == 0
                % "if the site above is also a"
                if isempty(find(sum([row_A col_A] == [g(site(1)) site(2)],2)==2,1)) == 0 % "if the site below is also a"
                    del_ind = find(sum([row_AA_v col_AA_v] == [g(site(1)) site(2)],2)==2,1);
                    row_AA_v(del_ind) = [];
                    col_AA_v(del_ind) = [];
                    rAA_v = rAA_v - kr;
                     
                     
                    row_AO_down = [row_AO_down; g(site(1))];
                    col_AO_down = [col_AO_down; site(2)];
                    rAO_down = rAO_down + kdiff;
                else
                    del_ind = find(sum([row_AO_up col_AO_up] == [site(1) site(2)],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                end
                
                % "if the site to the right is also a"
                if isempty(find(sum([row_A col_A] == [site(1) f(site(2))],2)==2,1)) == 0 % "if the site to the right is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) site(2)],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    
                    row_AO_left = [row_AO_left; site(1)];
                    col_AO_left = [col_AO_left; f(site(2))];
                    rAO_left = rAO_left + kdiff;
                else
                    del_ind = find(sum([row_AO_right col_AO_right] == [site(1) site(2)],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                end
                
                % If the site to the left is a-occupied, aa pair is annhiliated
                if isempty(find(sum([row_A col_A] == [site(1) g(site(2))],2)==2,1)) == 0 % "if the site to the left is also a"
                    del_ind = find(sum([row_AA_h col_AA_h] == [site(1) g(site(2))],2)==2,1);
                    row_AA_h(del_ind) = [];
                    col_AA_h(del_ind) = [];
                    rAA_h = rAA_h - kr;
                    
                    
                    row_AO_right = [row_AO_right; site(1)];
                    col_AO_right = [col_AO_right; g(site(2))];
                    rAO_right = rAO_right + kdiff;
                    
                else
                    del_ind = find(sum([row_AO_left col_AO_left] == [site(1) site(2)],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                end
            
                % If the site to the right of the diff. accept. is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [f(site(1)) f(site(2))],2)==2,1)) == 0 % "if the site to the right is a"
                    row_AA_h = [row_AA_h; f(site(1))];
                    col_AA_h = [col_AA_h; site(2)];
                    rAA_h = rAA_h + kr;
                    
                    del_ind = find(sum([row_AO_left col_AO_left] == [f(site(1)) f(site(2))],2)==2,1);
                    row_AO_left(del_ind) = [];
                    col_AO_left(del_ind) = [];
                    rAO_left = rAO_left - kdiff;
                else
                    row_AO_right = [row_AO_right; f(site(1))];
                    col_AO_right = [col_AO_right; site(2)];
                    rAO_right = rAO_right + kdiff;
                end
                
                % If the site to the left of the diff. accept. is a-occupied, aa pair is generated
                if isempty(find(sum([row_A col_A] == [f(site(1)) g(site(2))],2)==2,1)) == 0 % "if the site to the right is a"
                    row_AA_h = [row_AA_h; f(site(1))];
                    col_AA_h = [col_AA_h; g(site(2))];
                    rAA_h = rAA_h + kr;
                    
                    
                    del_ind = find(sum([row_AO_right col_AO_right] == [f(site(1)) g(site(2))],2)==2,1);
                    row_AO_right(del_ind) = [];
                    col_AO_right(del_ind) = [];
                    rAO_right = rAO_right - kdiff;
                else
                    row_AO_left = [row_AO_left; f(site(1))];
                    col_AO_left = [col_AO_left; site(2)];
                    rAO_left = rAO_left + kdiff;
                end
                
                % if the site below the diff. acceptor is also a
                if isempty(find(sum([row_A col_A] == [f(f(site(1))) site(2)],2)==2,1)) == 0 % "if the site below is a"
                    row_AA_v = [row_AA_v; f(site(1))];
                    col_AA_v = [col_AA_v; site(2)];
                    rAA_v = rAA_v + kr;
                    
                    
                    del_ind = find(sum([row_AO_up col_AO_up] == [f(f(site(1))) site(2)],2)==2,1);
                    row_AO_up(del_ind) = [];
                    col_AO_up(del_ind) = [];
                    rAO_up = rAO_up - kdiff;
                else
                    row_AO_down = [row_AO_down; f(site(1))];
                    col_AO_down = [col_AO_down; site(2)];
                    rAO_down = rAO_down + kdiff;
                end
            end
                
            % delete a and o
            del_ind = find(sum([row_A col_A] == [site(1) site(2)],2)==2,1);
            row_A(del_ind) = []; col_A(del_ind) = []; 
            del_ind = find(sum([row_O col_O] == [f(site(1)) site(2)],2)==2,1);
            row_O(del_ind) = []; col_O(del_ind) = []; 
            
            % add a and o
            row_O = [row_O; site(1)]; col_O = [col_O; site(2)]; % update O
            row_A = [row_A; f(site(1))]; col_A = [col_A; site(2)]; % update A
            
            % delete ao pair
            del_ind = find(sum([row_AO_down col_AO_down] == [site(1) site(2)],2)==2,1);
            row_AO_down(del_ind) = []; col_AO_down(del_ind) = [];
            rAO_down = rAO_down - kdiff;
            
            % make ao pair
            row_AO_up = [row_AO_up; f(site(1))]; col_AO_up = [col_AO_up; site(2)]; % update AO
            rAO_up = rAO_up + kdiff;
        end
    end
    %% Fractional Coverages
    thetaA = nA/d^2;
    thetaAA = nAA/d^2/2;
    thetaO = nO/d^2;
    thetaAO = thetaA - thetaAA;
    thetaOO = thetaO - thetaAO;
    muAA = thetaAA./thetaA./thetaA;
    %% Extract rates, reversibilities
    % Ensemble-specific -sorption rates
    raa_ads = TON_aa_ads./time; raa_des = TON_aa_des./time;
    rao_ads = TON_ao_ads./time; rao_des = TON_ao_des./time;
    STY = TON./time*2/d^2;
    
    % Ensemble-specific -sorption rates per 1st-order model (RSA)
    rao_des_RSA = -kAd*(thetaAA-thetaAO)*d^2;
    rao_ads_RSA = kA*(thetaOO-thetaAO)*d^2;
    
    raa_des_RSA = kAd*(thetaAA)*d^2;
    raa_ads_RSA = kA*(thetaAO)*d^2;
    
    % averaging post steady state
    em = round(length(thetaA)*0.25);
    
    % Reversibilities per kMC and 1st-order (RSA)
    zaa = mean(raa_des(end-em:end))./mean(raa_ads(end-em:end));
    zao = mean(rao_des(end-em:end))./mean(rao_ads(end-em:end));
    
    zaa_RSA = mean(raa_des_RSA(end-em:end))./mean(raa_ads_RSA(end-em:end));
    zao_RSA = mean(rao_des_RSA(end-em:end))./mean(rao_ads_RSA(end-em:end));
    
    % Putting repeats into arrays
    A(w) = mean(thetaA(end-em:end));
    AA(w) = mean(thetaAA(end-em:end));
    AO(w) = mean(thetaAO(end-em:end));
    MUAA(w) = mean(muAA(end-em:end));
    O(w) = mean(thetaO(end-em:end));
    R(w) = mean(STY(end-em:end));
    ZAA(w) = zaa;
    ZAO(w) = zao;
    
    % Make sure at steady state
%     figure
%     set(gca,'fontweight','bold','fontsize',11,'box','on');
%     set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
%     hold on
%     plot(time,STY,'k','linewidth',1.5)
%     plot(time(end-em:end),STY(end-em:end),'ko','linewidth',1.5)
%     xlabel('Poisson time / a.u.');
%     ylabel('TOF');
%     legend('kMC TOF','ss kMC TOF');
%     legend boxoff
end
%% Reaction rate
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(time,STY,'k','linewidth',1.5)
plot(time(end-em:end),STY(end-em:end),'ko','linewidth',1.5)
xlabel('Poisson time / a.u.');
ylabel('TOF');
legend('kMC TOF','ss kMC TOF');
legend boxoff
%% z adsorption for aa
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(time,raa_des_RSA./raa_ads_RSA,'k-','linewidth',0.5)
plot(time,raa_des./raa_ads,'r-','linewidth',1.5)
xlabel('Poisson time / a.u.');
ylabel('z_{ads,[{\itaa}]}');
legend('1^{st}-order eqn.','kMC');
legend boxoff
%% z adsorption for ao
figure
set(gca,'fontweight','bold','fontsize',11,'box','on');
set(gca,'linewidth', 2,'fontsize',12,'fontname','Arial')
hold on
plot(time,rao_des_RSA./rao_ads_RSA,'k-','linewidth',0.5)
plot(time,rao_des./rao_ads,'r-','linewidth',1.5)
xlabel('Poisson time / a.u.');
ylabel('z_{ads,[{\itao}]}');
legend('1^{st}-order eqn.','kMC');
legend boxoff
%% Extracting to command window and excel file
% M = [A;AA];
% writematrix(M,'M.xlsx')
% tock = [toc];
% [mean(A);std(A);mean(AA);std(AA);std(AA)/mean(AA)*100;mean(AA)/mean(A)^2;mean(ZAA);mean(ZAO)]
[mean(R);std(R)]
toc
