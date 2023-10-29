function [beta, pval, omega, tscore] = svcm(imgData, xMatrix)

lib_dir = '/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/OtherMethods/lib';
addpath(lib_dir)

ch = 1.1;
S = 10;
f1 = 1;
f2 = 1;

imgDataSize = size(imgData);
U = imgDataSize(1);
V = imgDataSize(2);
W = imgDataSize(3);
n = imgDataSize(4);
dimn = [U, V, W];
imgData = reshape(imgData, prod(imgDataSize(1:3)), imgDataSize(4));

[vxlNbySeqid, vxlDistSeq, vxlNberSeq] = findVoxelSequence(ch, S, dimn);
[mxBeta, mxCovb, mxR] = mlrWresd(xMatrix, imgData');
[smxR] = smoothResidual(mxR, dimn);
[mxAspBeta, mxAspCovb] = smoothParameter(vxlNbySeqid, vxlDistSeq, vxlNberSeq, mxR, smxR, mxBeta, xMatrix,S, f1, f2);

tscore = abs(mxAspBeta./mxAspCovb.^(0.5));
tdf = size(xMatrix,1) - size(xMatrix,2);
pval = 2 * (1 - tcdf(tscore, tdf));
% pSn = pval < 0.05;

p = size(mxAspBeta, 2);
beta = reshape(mxAspBeta, U, V, W, p);
pval = reshape(pval, U, V, W, p);
omega = reshape(smxR, U, V, W, n);
tscore = reshape(tscore, U, V, W, p);

end
