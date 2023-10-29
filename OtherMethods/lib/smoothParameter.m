function [fnlBeta test] = smoothParameter(vxlNbySeqid,vxlDistSeq, vxlNberSeq, mxR,smxR,mxBeta,xMatrix,hh,f1,f2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SVCM Project EXample1 June 30, 2011  LLK@CH %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input vxlNbySeqid
%   

sResidual = smxR';
inlBeta = mxBeta';
xDesign = xMatrix;


%% calculate the measurement error variance

mError = mxR - sResidual';
sigmaError = (std(mError)).^2;
sigmaError = sigmaError';


%% adaptively smooth parameters
ch = 1.1;
ss = 1:hh;
hSeq = ch.^ss;
hSeq = hSeq';

chiSeq = f1*chi2inv((0.80)./(ss(4:end)-2),1);
chiSeq = chiSeq';

% size(inlBeta) % N by p
% size(sResidual) % N by n
% size(hSeq) % nh by 1
% size(vxlNberSeq) % N by 1
% size(vxlNbySeqid) % N by m
% size(vxlDistSeq) % N by m
% size(xDesign) % n by p
% size(chiSeq) % nh-3 by 1
% size(sigmaError) % N by 1

n = size(xDesign,1);
Cn = f2*log(n)*chi2inv(0.95,1);

[fnlBeta OptHs test] = aspV2(inlBeta,sResidual,hSeq,vxlNberSeq,vxlNbySeqid,vxlDistSeq,xDesign,chiSeq,sigmaError,Cn);

end