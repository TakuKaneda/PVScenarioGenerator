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
sun_range = min_hour:max_hour;

toc;
disp('sun rise & set done');
disp(['min sunrize:',num2str(min_hour) ,' max sunset: ',num2str(max_hour)]);
% save sunrize & sunset data
sunname = ['pdfdata/', datafilename,'_',num2str(N),'_sundata.mat'];
save(sunname,'sunx1','sunx2','fsun','sun_range','H');
clear 'sunx1' 'sunx2' 'fsun'
%% Step 2 : Generate the joint probability
disp('start generation of joint pdfs');
% n_rep = 4; % number of repeat of the algorithm
% for i = 1:n_rep
%     xx = [floor(length(v)/n_rep*(i-1)+1), floor(length(v)/n_rep*i)]
% end
for i = 1:n_rep
    s_time = sun_range(floor(length(sun_range)/n_rep*(i-1)+1));
    e_time = sun_range(floor(length(sun_range)/n_rep*i));
    pv_jointpb =zeros(N,N,e_time-s_time+1);
    p_range1=zeros(N,e_time-s_time+1);
    p_range2=zeros(N,e_time-s_time+1);
    for t = 1:e_time-s_time+1
        % 2.1 : take hi2 : bandwidth of the murlivariate kernel estimatiion
        %       and   x  : range of PV power for each time
        [~,x,hi2] = ksdensity(pv_data(:,(s_time-1)+t-1:(s_time-1)+t));
        x1 = reshape(x(:,1),[30 30]);   
        x2 = reshape(x(:,2),[30 30]);
        x1 = x1(1,:);
        x2 = x2(:,1);

        % 2.2 : estimate the joint PDF of i and i-1 by KDF
        p_range1(:,t) = linspace(x1(1),x1(end),N); % pv range of time t-1
        p_range2(:,t) = linspace(x2(1),x2(end),N); % pv range of time t

        hi2 = hi2*0.5; % HALF
        pv_jointpb(:,:,t) = joint_pdf(hi2,n_days,p_range1(:,t),p_range2(:,t),pv_data,(s_time-1)+t);
        toc;
        disp(['time ', num2str((s_time-1)+t) ,' done']);
    end
    disp([num2str(i),'/',num2str(n_rep),' has finished'])
    pdfdata_name = ['pdfdata/', datafilename,'_',num2str(N),'_',num2str(i),'_jointpdf.mat'];
    save(pdfdata_name,'s_time','e_time','p_range1','p_range2','pv_jointpb');
end
