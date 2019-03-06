%% 1 - Building_Prob_run;
% tic();
N = 100; % mesh size, recommended to be at least 1000.
% It takes time to compute but it must only be computed once for a given
% data set to catch the correlation between adjacent hours.

% datafilename = '15m_157_9.99_winter';
datafilename = 'PGE-SASH-4101_winter';
pv_data =  csvread(['../preprocessing/',datafilename,'.csv']);

% permutatoin and create train & test sets
train_rate = .75; % < 1
pv_data = pv_data(randperm(length(pv_data)),:);
perm_idx = randperm(length(pv_data));
train_idx = 1:floor(length(pv_data)*train_rate);
test_idx = train_idx(end)+1:length(pv_data);

% Computation of joint PDFs (PV power of time t and t-1 / Sun rise and sun set)
n_rep = 20; % NUMBER OF DATA SPLIT
disp('Starting to build the probability density functions');
pv_probability_decomposed(pv_data(train_idx,:),N,n_rep,datafilename);

%% 2 - Sampling_run;
% Variables
n_scenario = 500; % # of scenario
% The computation is quick, can generate any number of scenarios from the
% PDF computed at step 1.

% Sampling
disp('Starting to generate the scenarios');
[pv_scenario,tss,tes]=pv_sampling_decomposed(n_scenario,N,n_rep,datafilename);
filename = ['../pv_simulation/',datafilename,'_',num2str(n_scenario),'sim_',num2str(N),'meshpdf.csv'];
% csvwrite(filename,pv_scenario)
%% 3 - Visualization
timesteps = linspace(0,23.75,size(pv_scenario,2));
figure()
subplot(2,2,1);
plot(timesteps,pv_data(test_idx,:));
xlabel('time (h)')
ylabel('power (kW)')
title('test data');

subplot(2,2,2);
plot(timesteps,pv_scenario);
xlabel('time (h)')
ylabel('power (kW)')
title('simulated data');

subplot(2,2,3);
hold on
plot(timesteps,mean(pv_data(test_idx,:)))
plot(timesteps,mean(pv_scenario));
% errorbar(timesteps,mean(pv_data),std(pv_data));
% errorbar(timesteps,mean(pv_scenario),std(pv_scenario));
title('mean production');
xlabel('time (h)')
ylabel('power (kW)')
legend('test data','simulation');

subplot(2,2,4);
[f,xi] = ksdensity(sum(pv_data(test_idx,:),2)./4);
[f2,xi2] = ksdensity(sum(pv_scenario,2)./4);
plot(xi,f);
hold on
plot(xi2,f2);
title('density of production per day')
xlabel('energy (kWh)')
legend('test data','simulation');
% For density of power at a certain timestep
% hour = 13;
% [f,xi] = ksdensity(pv_data(:,hour*4));
% [f2,xi2] = ksdensity(pv_scenario(:,hour*4));
% plot(xi,f);
% hold on
% plot(xi2,f2);
% title(['density of production at a hour ',num2str(hour)]);
% xlabel('power (kW)')

%% 4 - Error analysis
% compute corelation coefficient
R = corrcoef(pv_scenario);
Rreal = corrcoef(pv_data(test_idx,:));
% adjacent relation
R_adj = diag(R,1);
Rreal_adj = diag(Rreal,1);
error = (Rreal-R)./Rreal*100; % error matrix
error_adj = diag(error,1); % adjacent error