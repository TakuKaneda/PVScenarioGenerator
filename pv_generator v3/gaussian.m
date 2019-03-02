function y = gaussian(x)
% Authors : Bruno Losseau <bruno.losseau@student.uclouvain.be>
%         : Taku Kaneda <taku.kaneda@student.uclouvain.be>
% Course  : LINMA2360 Project in mathematical engineering
% Date    : 19th Nov. 2016

% this function compute 'standard' Gaussian distribution (or Normal distribution)
% 'standard' means that mean = 0, variance = 1
% input  - x : random variable
% output - y : probability of taking x
    y = 1./sqrt(2.*pi).*exp(-1/2.*x.^2);
end