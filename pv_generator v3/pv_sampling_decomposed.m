function [pv_scenario,tss,tes]=pv_sampling_decomposed(n_scenario,N,n_rep,pv_capacity,datafilename)
tic;
%% load sun data
sun_load = load(['pdfdata/', datafilename,'_',num2str(N),'_sundata.mat']);
fsun = sun_load.fsun; H = sun_load.H;
sunx1 = sun_load.sunx1; sunx2 = sun_load.sunx2;
sun_range = sun_load.sun_range;
clear sun_load

pv_scenario = zeros(n_scenario,H);
% generate sun rize/set scenarios 
maxSun = max(max(fsun));
tss = zeros(n_scenario,1);
tes = zeros(n_scenario,1);
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
end

clear fsun sunx1 sunx2 sun_range
%% generate PV scenarios
tic
for i = 1:n_rep
    % load the i-th joint pdf data
    jpdf_load = load(['pdfdata/', datafilename,'_',num2str(N),'_',num2str(i),'_jointpdf.mat']);
    s_time = jpdf_load.s_time ; e_time = jpdf_load.e_time ;
    pv_jointpb = jpdf_load.pv_jointpb;
    p_range1 = jpdf_load.p_range1; p_range2 = jpdf_load.p_range2;
    clear jpdf_load
    
    for n = 1:n_scenario
        % 2.3 : generate one random sample of PV power at i by rejection sampling method
        stime_n = s_time; etime_n = e_time;
        if i == 1
            stime_n = tss(n); 
        elseif i == n_rep
            etime_n = tes(n);
        end
        if etime_n > sun_range(end)  % if sun set time is greater than the data -> to be max val of data
            etime_n = sun_range(end);
        end
        for t = stime_n:etime_n % real time loop
            tt = t-stime_n+1;
            [~,pre_pv_index] = min(abs(pv_scenario(n,t-1)-p_range1(:,tt)));    % take index of PV power ar i-1
            cumul_proba = cumulProba(pv_jointpb(:,pre_pv_index(1),tt),p_range2(:,tt));    % make conditional CDF in terms of the value of PV power at i-1
            random = rand; % random value in [0,1]
            [~,pv_index] = min(abs(cumul_proba-random));

           if p_range2(pv_index,tt) < 0 % judge : lower bound (=0,nonnegative)
               pv_scenario(n,t) = 0; 
           else
               if p_range2(pv_index,tt) > pv_capacity
                   pv_scenario(n,t) = pv_capacity; % judge : upper bound (=pv_capacity)
               else    
                   pv_scenario(n,t) = p_range2(pv_index,tt);
               end
           end
        end
    end
    disp([num2str(floor(i/n_rep*100)),'% is done.']);
toc
end
end
