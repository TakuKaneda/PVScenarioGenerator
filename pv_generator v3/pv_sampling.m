function [pv_scenario,tss,tes]=pv_sampling(n_scenario,pv_capacity,sunx1,sunx2,fsun,p_range1,p_range2,pv_jointpb)
tic;
n_hours = length(p_range1(1,:));
N = length(p_range1(:,1));
pv_scenario = zeros(n_scenario,n_hours);
maxSun = max(max(fsun));
tss = zeros(n_scenario,1);
tes = zeros(n_scenario,1);
tic
for n = 1:n_scenario
    % 1.3 : take one random sample of (tss,tes) by rejection sampling method
    bin = 0;
    while ~bin
        uniform_random_pv = unifrnd(0,maxSun);
        ts =  randi(N);
        te =  randi(N);
        if uniform_random_pv < fsun(ts,te)
            bin=1;
        end
    end
    tss(n) = round(sunx1(ts));
    tes(n) = round(sunx2(te));
    
    % 2.3 : generate one random sample of PV power at i by rejection sampling method
    for i = tss(n):tes(n)
        [~,pre_pv_index] = min(abs(pv_scenario(n,i-1)-p_range1(:,i)));    % take index of PV power ar i-1
        cumul_proba = cumulProba(pv_jointpb(:,pre_pv_index(1),i),p_range2(:,i));    % make conditional CDF in terms of the value of PV power at i-1
        random = rand; % random value in [0,1]
        [~,pv_index] = min(abs(cumul_proba-random));
        
       if p_range2(pv_index,i) < 0 % judge : lower bound (=0,nonnegative)
           pv_scenario(n,i) = 0; 
       else
           if p_range2(pv_index,i) > pv_capacity
               pv_scenario(n,i) = pv_capacity; % judge : upper bound (=pv_capacity)
           else    
               pv_scenario(n,i) = p_range2(pv_index,i);
           end
       end
    end
    if mod(n,1000)==999
        disp(['scenario ', num2str(n+1), ' done']);
        toc;
    end
end
toc
end
