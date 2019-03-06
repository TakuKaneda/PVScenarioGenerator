function [pv_scenario,tss,tes]=pv_sampling_decomposed(n_scenario,N,n_rep,datafilename)
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
        % exit if the value pair is under the pdf and falls in the range of
        % sunrise/set time of the data
        if uniform_random_pv < fsun(ts,te) && ...
                round(sunx1(ts)) >= sun_range(1) && ...
                round(sunx2(te)) <= sun_range(end)
            tss(n) = round(sunx1(ts));
            tes(n) = round(sunx2(te));
            bin=1;
        end
    end
%     tss(n) = round(sunx1(ts));
%     tes(n) = round(sunx2(te));
%     
%     % checking the range
%     if tss(n) < sun_range(1) % if sunrise time is smaller than the data -> to be min val of data
%         tss(n) = sun_range(1);
%     end
%     if tes(n) > sun_range(end) % if sunset time is greater than the data -> to be max val of data
%         tes(n) = sun_range(end);
%     end
end

clear fsun sunx1 sunx2
%% generate PV scenarios
tic
for i = 1:n_rep
    % load the i-th joint pdf data
    jpdf_load = load(['pdfdata/', datafilename,'_',num2str(N),'_',num2str(i),'_jointpdf.mat']);
    s_time = jpdf_load.s_time ; e_time = jpdf_load.e_time ;
    pv_jointpb = jpdf_load.pv_jointpb;
    p_range1 = jpdf_load.p_range1; p_range2 = jpdf_load.p_range2;
    pv_range = jpdf_load.pv_range;
    clear jpdf_load
    
    idx = find(tss<=e_time & s_time<=tes); % indices to be simulated
    for m = 1:length(idx) 
        n = idx(m);
        % 2.3 : generate one random sample of PV power by rejection sampling method
        n_s_time = s_time; n_e_time = e_time;
        if s_time < tss(n)
            n_s_time = tss(n); 
        elseif tes(n) < e_time
            n_e_time = tes(n);
        end
        for t = n_s_time:n_e_time % t: real time loop
            tt = t-n_s_time+1; % tt: temporal time index
            [~,pre_pv_index] = min(abs(pv_scenario(n,t-1)-p_range1(:,tt)));    % take index of PV power at t-1
            cumul_proba = cumulProba(pv_jointpb(:,pre_pv_index(1),tt),p_range2(:,tt));    % make conditional CDF in terms of the value of PV power at t-1
            
            % sampling of the current PV output
            bin = 0;
            while ~bin
                [~,pv_index] = min(abs(cumul_proba-rand)); % index of the p_range2
                % exit if the power is positive & falls within the data range 
                if p_range2(pv_index,tt) > 0 && ...
                        pv_range(tt,1) <= p_range2(pv_index,tt) && ...
                        p_range2(pv_index,tt) <= pv_range(tt,2)
                    pv_scenario(n,t) = p_range2(pv_index,tt); 
                    bin = 1;
                end
            end
            
        end
    end
    disp([num2str(floor(i/n_rep*100)),'% done.']);
toc
end
end
