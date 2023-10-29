function [img xMatrix bet et] = svcm_gendata(n, noi, seed)

% This script can only generate
% img dim: 64 x 64 x 8
% num covariates: 3

rng(seed);

%lib_dir = 'lib';
addpath(lib_dir)

pattern = 1; % one of 1:7
pr = 1; % controls beta (main effect) scales
nf1 = 1; % controls eta (individual effect) scales
ef = 1; % 1 for xi=norm(0,1), 2 for x=chisq(3)-3, where xi is the n x 3 matrix to be eigendcomposed to generate the coefficients of individual effects

xMatrix = ones(n,3);
xMatrix(:,2) = binornd(1,0.5,n,1);
xMatrix(:,3) = unifrnd(1,2,n,1);
xMatrix(:,2) = (xMatrix(:,2)-mean(xMatrix(:,2)))/std(xMatrix(:,2));
xMatrix(:,3) = (xMatrix(:,3)-mean(xMatrix(:,3)))/std(xMatrix(:,3));

[beta roiIds] = generateDesign(pattern, pr); % can only generate 64 x 64 x 3 designs
[bet et img] = generateImageData(beta, xMatrix, nf1, ef, noi, seed); % stack betas 8 times
end


