function biv = joint_pdf(hi2,n_days,p_range1,p_range2,pv_data,time)
% Authors : Bruno Losseau <bruno.losseau@student.uclouvain.be>
%         : Taku Kaneda <taku.kaneda@student.uclouvain.be>
% Course  : LINMA2360 Project in mathematical engineering
% Date    : 19th Nov. 2016

% This function compute the Multivariate kernel density estimation of PV power at i and i-1
% Gaussian function is used as the kernel function
% input - hi2        : bandwidth of the murlivariate kernel estimatiion
%       - n_days     : number of days in data
%       - p_range1   : range of PV power at i-1
%       - p_range2   : range of PV power at i
%       - pv_24hours : real data
%       - time       : number of hour (=i)
% output- biv       : value of joint probability in matrix form

% tic
N = length(p_range1);
biv = zeros(N,N);

for k = 1:N
    for l = 1:N
        i = 1:n_days;
%         y = 1./n_days./hi2(2)./hi2(1).*gaussian((p_range2(l)-pv_24hours(i,time))./hi2(2))...
%                  .*gaussian((p_range1(k)-pv_24hours(i,time-1))./hi2(1));
        y = 1./n_days./hi2(2)./hi2(1).*normpdf((p_range2(l)-pv_data(i,time))./hi2(2))...
                 .*normpdf((p_range1(k)-pv_data(i,time-1))./hi2(1));
        biv(l,k) = sum(y);
    end
end
% toc

end