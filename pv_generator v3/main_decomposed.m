%% 1 - Building_Prob_run;
% tic();
N = 4000; % mesh size, recommended to be at least 1000.
% It takes time to compute but it must only be computed once for a given
% data set to catch the correlation between adjacent hours.
datafilename = '15m_157_9.99_winter';
pv_data =  csvread(['../preprocessing/',datafilename,'.csv']);
% Computation of joint PDFs (PV power of time t and t-1 / Sun rise and sun set)
n_rep = 20; % NUMBER OF DATA SPLIT
disp('Starting to build the probability density functions');
% % pv_probability_decomposed(pv_data,N,n_rep,datafilename);

%% 2 - Sampling_run;
% Variables
n_scenario = 100; % # of scenario
% The computation is quick, can generate any number of scenarios from the
% PDF computed at step 1.
pv_capacity = 9.99;  % capacity of the pv generator
% H = size(pv_data,2); % horizon
% sundata_name = ['pdfdata/', datafilename,'_',num2str(N),'_sundata.mat'];
% pdfdata_name = 
% Sampling
disp('Starting to generate the scenarios');
% [pv_scenario,tss,tes]=pv_sampling(n_scenario,pv_capacity,sunx1,sunx2,fsun,p_range1,p_range2,pv_jointpb);
% % [pv_scenario,tss,tes]=pv_sampling_decomposed(n_scenario,N,n_rep,pv_capacity,datafilename);
% filename = ['../pv_simulation/',datafilename,'_',...
%             num2str(n_scenario),'sim_',num2str(N),'meshpdf.csv'];
% csvwrite(filename,pv_scenario)
%% 3 - Error analysis
timesteps = 1:size(pv_scenario,2);
plot(timesteps,pv_scenario);
% compute corelation coefficient
R = corrcoef(pv_scenario);
Rreal = corrcoef(pv_data);
% adjacent relation
R_adj = diag(R,1);
Rreal_adj = diag(Rreal,1);
error = (Rreal-R)./Rreal*100; % error matrix
error_adj = diag(error,1); % adjacent error