function [sunx1,sunx2,fsun,p_range1,p_range2,pv_jointpb]=pv_probability(pv_data,N)
%% define the parameters 
n_days =length(pv_data(:,1)); % number of days in data
n_hours = length(pv_data(1,:)); % number of hours per day in data
%% Step 1 : Generate random samples 'tss'(=time of sunrise) and 'tes'(=time of sunset)
% 1.1 : find the time of sunrise and sunset in data(pv_24hours)
tic;
disp('start sun rise & set');
first_sun = zeros(n_days,1);  
last_sun = zeros(n_days,1); 
for t=1:n_days
    first_sun(t) = find(pv_data(t,:),1);
    last_sun(t) = find(pv_data(t,:),1,'last');
end
sun_hours = [first_sun(:,1) last_sun(:,1)];

min_hour = min(first_sun(:,1));
max_hour = max(last_sun(:,1));

% 1.2 : estimate bandwidht of KDF and ranges(domain)
[~,xsun,hi1] = ksdensity(sun_hours);
sunrange1=reshape(xsun(:,1),[30 30]);
sunrange2=reshape(xsun(:,2),[30 30]);
sunx1 = sunrange1(1,:);
sunx2 = sunrange2(:,1);
 
% Compute the PDF
sunx1 = linspace(sunx1(1),sunx1(end),N);
sunx2 = linspace(sunx2(1),sunx2(end),N);

hi1 = hi1.*0.5; %HALF
fsun = joint_pdf(hi1,n_days,sunx1,sunx2,sun_hours,2);

toc;
disp('sun rise & set done');
%% Step 2 : Generate the joint probability
disp('start generation of joint pdfs');
pv_jointpb =zeros(N,N,n_hours);
p_range1=zeros(N,n_hours);
p_range2=zeros(N,n_hours);

for t = min_hour:max_hour    
    % 2.1 : take hi2 : bandwidth of the murlivariate kernel estimatiion
    %       and   x  : range of PV power for each time
    [~,x,hi2] = ksdensity(pv_data(:,t-1:t));
    x1 = reshape(x(:,1),[30 30]);   
    x2 = reshape(x(:,2),[30 30]);
    x1 = x1(1,:);
    x2 = x2(:,1);

    % 2.2 : estimate the joint PDF of i and i-1 by KDF
    p_range1(:,t) = linspace(x1(1),x1(end),N); % pv range of time i-1
    p_range2(:,t) = linspace(x2(1),x2(end),N); % pv range of time i
    
    hi2 = hi2*0.5; % HALF
    pv_jointpb(:,:,t) = joint_pdf(hi2,n_days,p_range1(:,t),p_range2(:,t),pv_data,t);
    toc;
    disp(['time ', num2str(t) ,' done']);
end

end
