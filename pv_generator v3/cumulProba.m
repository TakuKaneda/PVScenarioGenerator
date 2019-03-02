function y = cumulProba(f,range)
% Authors : Bruno Losseau <bruno.losseau@student.uclouvain.be>
%         : Taku Kaneda <taku.kaneda@student.uclouvain.be>
% Course  : LINMA2360 Project in mathematical engineering
% Date    : 19th Nov. 2016

% This function computes cumulative distribution function (CDF) by using
% probability density function (PDF) and the value of random variable
% input  - f     : PDF
%        - range : value of random variable
% output - y     : CDF

% define the parameters
N = length(range);
proba =f/sum(f); %normalise the value of PDF
cumul_proba = [proba(1) zeros(1,N-1)];

% compute CDF
for j=2:N
    cumul_proba(j) = cumul_proba(j-1)+proba(j);
end
y = cumul_proba;
end