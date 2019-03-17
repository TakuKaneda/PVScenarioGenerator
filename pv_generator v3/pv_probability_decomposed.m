function pv_probability_decomposed(pv_data,N,n_rep,datafilename)
%% define the parameters 
n_days =length(pv_data(:,1)); % number of days in data
H = length(pv_data(1,:)); % number of time steps per day in data
%% Step 1 : Generate random samples 'tss'(=time of sunrise) and 'tes'(=time of sunset)
% 1.1 : find the time of sunrise and sunset in data(pv_24hours)
tic;
disp('start sun rise & set');
first_sun = zeros(n_days,1);  
last_sun = zeros(n_days,1); 
for tt=1:n_days
    first_sun(tt) = find(pv_data(tt,:),1);
    last_sun(tt) = find(pv_data(tt,:),1,'last');
end
sun_hours = [first_sun(:,1) last_sun(:,1)];

min_hours = min(sun_hours);
max_hours = max(sun_hours);

% 1.2 : estimate bandwidht of KDF and ranges(domain)
[~,~,hi1] = ksdensity(sun_hours);

% % Compute the PDF
% sunx1 = linspace(xsun(1,1),xsun(end,1),N);
% sunx2 = linspace(xsun(1,2),xsun(end,2),N);
sunx1 = linspace(min_hours(1),max_hours(1),N);
sunx2 = linspace(min_hours(2),max_hours(2),N);

half_hi1 = hi1.*0.4; %scale
fsun = joint_pdf(half_hi1,n_days,sunx1,sunx2,sun_hours,2);
sun_range = min_hours(1):max_hours(2);

toc;
disp('sun rise & set done');
disp(['min sunrize:',num2str(min_hours(1)) ,' max sunset: ',num2str(max_hours(2))]);
% save sunrize & sunset data
sunname = ['pdfdata/', datafilename,'_',num2str(N),'_sundata.mat'];
save(sunname,'sunx1','sunx2','fsun','sun_range','H');
clear 'sunx1' 'sunx2' 'fsun'
%% Step 2 : Generate the joint probability
disp('start generation of joint pdfs');
for i = 1:n_rep
    s_time = sun_range(floor(length(sun_range)/n_rep*(i-1)+1));
    e_time = sun_range(floor(length(sun_range)/n_rep*i));
    pv_jointpb =zeros(N,N,e_time-s_time+1);
    p_range1=zeros(N,e_time-s_time+1);
    p_range2=zeros(N,e_time-s_time+1);
    pv_range = zeros(e_time-s_time+1,2);
    for tt = 1:e_time-s_time+1
        % 2.1 : take hi2 : bandwidth of the murlivariate kernel estimatiion
        %       and   x  : range of PV power for each time
        [~,~,hi2] = ksdensity(pv_data(:,(s_time-1)+tt-1:(s_time-1)+tt));
%         x1 = reshape(x(:,1),[30 30]);   
%         x2 = reshape(x(:,2),[30 30]);
%         x1 = x1(1,:);
%         x2 = x2(:,1);

        % 2.2 : estimate the joint PDF of i and i-1 by KDF
        x1_start = min(pv_data(:,(s_time-1)+tt-1)); 
        x1_end = max(pv_data(:,(s_time-1)+tt-1));  
        x2_start = min(pv_data(:,(s_time-1)+tt)); 
        x2_end = max(pv_data(:,(s_time-1)+tt)); 
%         x1_start = x(1,1); x1_end = x(end,1);
%         x2_start = x(1,2); x2_end = x(end,2);
        
        p_range1(:,tt) = linspace(x1_start,x1_end,N); % pv range of time t-1
        p_range2(:,tt) = linspace(x2_start,x2_end,N); % pv range of time t

        half_hi2 = hi2*0.5; % HALF
        pv_jointpb(:,:,tt) = joint_pdf(half_hi2,n_days,p_range1(:,tt),p_range2(:,tt),pv_data,(s_time-1)+tt);
        
        pv_range(tt,1) = min(pv_data(:,(s_time-1)+tt));
        pv_range(tt,2) = max(pv_data(:,(s_time-1)+tt));
        toc;
        disp(['time ', num2str((s_time-1)+tt) ,' done']);
    end
    disp([num2str(i),'/',num2str(n_rep),' has finished'])
    pdfdata_name = ['pdfdata/', datafilename,'_',num2str(N),'_',num2str(i),'_jointpdf.mat'];
    save(pdfdata_name,'s_time','e_time','p_range1','p_range2','pv_jointpb','pv_range');
end
