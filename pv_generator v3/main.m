%% 1 - Building_Prob_run;
% tic();
N = 200; % mesh size, recommended to be at least 1000.
% It takes time to compute but it must only be computed once for a given
% data set to catch the correlation between adjacent hours.
datafilename = '157_9.99_winter';
pv_data =  csvread(strcat('../preprocessing/',datafilename,'.csv'));
% Computation of joint PDFs (PV power of time t and t-1 / Sun rise and sun set)
disp('Starting to build the probability density functions');
[sunx1,sunx2,fsun,p_range1,p_range2,pv_jointpb]=pv_probability(pv_data,N);
% toc();
pdfname = ['../pv_simulation/',datafilename,'_',num2str(N),'_jointpdf.mat'];
save(pdfname,'pv_jointpb');
%% 2 - Sampling_run;
% Variables
n_scenario = 1000; % # of scenario
% The computation is quick, can generate any number of scenarios from the
% PDF computed at step 1.
pv_capacity = 9.99;  % capacity of the pv generator
H = size(pv_data,2); % horizon

% Sampling
disp('Starting to generate the scenarios');
[pv_scenario,tss,tes]=pv_sampling(n_scenario,pv_capacity,sunx1,sunx2,fsun,p_range1,p_range2,pv_jointpb);

filename = ['../pv_simulation/',datafilename,'_',...
            num2str(n_scenario),'sim_',num2str(N),'meshpdf.csv'];
csvwrite(filename,pv_scenario)
%% 3 - Error analysis
% hours = linspace(1,H,H);
% plot(hours,pv_scenario);
% % compute corelation coefficient
% R = corrcoef(pv_scenario);
% Rreal = corrcoef(pv_data);
% % adjacent relation
% R_adj = diag(R,1);
% Rreal_adj = diag(Rreal,1);
% error = (Rreal-R)./Rreal*100; % error matrix
% error_adj = diag(error,1); % adjacent error