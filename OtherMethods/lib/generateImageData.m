function [bet, et, img] = generateImageData(beta, xMatrix, nf1, ef, noi, seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SVCM Project EXample1 June 29, 2011  LLK@CH %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a 64 by 64 by 8 3d image with n = 75

rng(seed);

n = size(xMatrix,1);

%% load image data and covariates


%% generate response image %%%
imgMatrixT = zeros(n,64,64);
for i=1:n
    imgMatrixT(i,:,:) = xMatrix(i,1)*beta(:,:,1) + xMatrix(i,2)*beta(:,:,2) + xMatrix(i,3)*beta(:,:,3);
end
%% generate response image %%%

%% generate response slice with error noises

phi1f = 'sqrt(1/4)*sin(2*pi*d1./64)';
phi2f = 'sqrt(1/4)*cos(2*pi*d2./64)';
phi3f = 'sqrt(1/2.625)*(9/8-d3./4)';

d1 = 1:64;
d2 = 1:64;
d3 = 1:8;

phi1 = eval(phi1f);
phi2 = eval(phi2f);
phi3 = eval(phi3f);

lambda1 = sqrt(.3*nf1);
lambda2 = sqrt(.15*nf1);
lambda3 = sqrt(.05*nf1);

sigma = sqrt(1);

imgData = zeros(64*64*8,n);

if ef==1
   xiMatrix = normrnd(0,1,n,3);
elseif ef==2
   xiMatrix = normrnd(0,1,n,3).^2 + normrnd(0,1,n,3).^2 + normrnd(0,1,n,3).^2 -3;
end

[VV DD EE] = pca(xiMatrix);
for i=1:3
    xiMatrix(:,i) = DD(:,i)/sqrt(EE(i));
end

if noi == 1
    noise = normrnd(0, sigma, 64, 64, 8, n);
elseif noi == 2
    noise = normrnd(0, sigma, 64, 64, 8, n).^2 + normrnd(0, sigma, 64, 64, 8, n).^2 + normrnd(0, sigma, 64, 64, 8, n).^2 - 3;
end

ph1 = repmat(phi1', 1, 64, 8) * lambda1;
ph2 = repmat(phi2, 64, 1, 8) * lambda2;
ph3 = permute(repmat(phi3', 1, 64, 64), [3 2 1 4]) * lambda3;
phi = [reshape(ph1, [], 1) reshape(ph2, [], 1) reshape(ph3, [], 1)];
et = phi * xiMatrix';
et = reshape(et, 64, 64, 8, n);
bet = permute(repmat(permute(beta, [3 1 2]), [1 1 1 8]), [2 3 4 1]);
img = reshape(bet, [], 3) * xMatrix' + reshape(et, [], n) + reshape(noise, [], n);
img = reshape(img, 64, 64, 8, n);

end
