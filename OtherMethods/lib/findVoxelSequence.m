function [VSeqId VSeqDist InVSeq VSeq] = findVoxelSequence(ch, S, dimn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SVCM Project EXample1 June 30, 2011  LLK@CH %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% revise @ June 15 2011
% treat the distance between voxel as 1
% redesign radius sequence 
% revised @ June 20, 2011
% sort voxel sequence in distance order
% changed the distance square to distance
% revised @ June 21, 2011
% applied WM only mask and add finding sequecne order in nozero id
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% revised to adapt to example 1
% VSeq - the oder of nearby voxel sequence
% VSeqDist - the distance of nearby voxel sequence
% InVSeq - maximum number of nearby vxoel
% VSeqId - Id of nearvy voxel

%% Given ch and S, find the nearby voxel sequence

% clear all;
% 
% homeDir1 = 'H:';
% homeDir2 = '/home/users/linglong';
% 
% homeDir = homeDir2;
% 
% dataPath = sprintf('%s/PDF/2012/Project/SVCM/2012_01_30_Example1/2012_03_07_Normal75/data',homeDir);
% 
% aa = sprintf('%s/vxlSeq.mat',dataPath);
% bb = sprintf('%s/vxlSeqDist.mat',dataPath); 
% cc = sprintf('%s/noVxlSeq.mat',dataPath);
% dd = sprintf('%s/vxlSeqId.mat',dataPath);

%% calculate the nearby voxel sequence
% dimn=[64 64 8];
% ch = 1.1;
% S = 10;
mask = ones(dimn);
maskVector = mask(:);
nonzeroId = find(maskVector~=0);
nId = numel(nonzeroId); 

[yy xx zz] = meshgrid(1:dimn(2),1:dimn(1),1:dimn(3));
xx = xx(nonzeroId); % xx row index of nonzero ID
yy = yy(nonzeroId); % yy column index of nonzero ID
zz = zz(nonzeroId); % zz depth index of nonzero ID

Lx = min(xx);
Ux = max(xx);
Ly = min(yy);
Uy = max(yy);
Lz = min(zz);
Uz = max(zz);

hS = ch^S;
hS2 = hS^2;
hSx = floor(hS -0.00001);
hSy = floor(hS -0.00001);
hSz = floor(hS -0.00001);

tic;
MnVSeq = 50;
VSeq = zeros(nId,MnVSeq);
VSeqDist = zeros(nId,MnVSeq);
InVSeq = zeros(nId,1);

for Nidi=1:nId
    TVSeq = [];
    TVSeqDist = [];
    
    tx = xx(Nidi);
    ty = yy(Nidi);
    tz = zz(Nidi);
      
    Ltx = max(tx-hSx,Lx);
    Utx = min(tx+hSx,Ux);
    Lty = max(ty-hSy,Ly);
    Uty = min(ty+hSy,Uy);
    Ltz = max(tz-hSz,Lz);
    Utz = min(tz+hSz,Uz);
    
    Ntx = Utx-Ltx+1;
    Nty = Uty-Lty+1;
    Ntz = Utz-Ltz+1;
    
    [Newty Newtx Newtz] = meshgrid(Lty:Uty,Ltx:Utx,Ltz:Utz);
    
    Newtx = Newtx(:);
    Newty = Newty(:);
    Newtz = Newtz(:);

    tdist = (tx-Newtx).^2+(ty-Newty).^2+(tz-Newtz).^2;
    tdistid = find(tdist<hS2);
    tdist = tdist(tdistid);
    SeqID = Newtx(tdistid) + (Newty(tdistid)-1)*dimn(1) + (Newtz(tdistid)-1)*dimn(1)*dimn(2);
    tMask = maskVector(SeqID);
     
     for ii=1:numel(tdistid)
         if (tMask(ii)~=0)
             TVSeq = [TVSeq SeqID(ii)];
             TVSeqDist = [TVSeqDist tdist(ii)];
         end             
     end 
     
     TnVSeq = numel(TVSeq);
     if (MnVSeq < TnVSeq) 
         VSeq = [VSeq zeros(nId,(TnVSeq-MnVSeq))];
         VSeqDist = [VSeqDist zeros(nId,(TnVSeq-MnVSeq))];
         MnVSeq = TnVSeq;
     end
             
     VSeq(Nidi,1:TnVSeq) = TVSeq;  
     VSeqDist(Nidi,1:TnVSeq) = TVSeqDist; 
     InVSeq(Nidi) = TnVSeq;
end

VSeq = VSeq(:,1:MnVSeq);
VSeqDist = VSeqDist(:,1:MnVSeq);

%% sort voxel sequence according to distance
for iId = 1:nId
    [Btemp, IXtemp] = sort(VSeqDist(iId,1:InVSeq(iId)));
    VSeqDist(iId,1:InVSeq(iId)) = Btemp;
    VSeq(iId,1:InVSeq(iId)) = VSeq(iId,IXtemp);
end
%% sort voxel sequence according to distance

%% transfer voxel sequence to be indexed by id
VSeqId = zeros(size(VSeq));
for iId = 1:nId
    subVSeq = VSeq(iId,1:InVSeq(iId));
    [tf, loc] = ismember(subVSeq,nonzeroId);
    VSeqId(iId,1:InVSeq(iId)) = loc;
end
%% transfer voxel sequence to be indexed by id

VSeqDist = VSeqDist.^0.5;

% save(aa,'VSeq','-mat')
% save(bb,'VSeqDist','-mat')
% save(cc,'InVSeq','-mat')
% save(dd,'VSeqId','-mat')
% 
% tt = toc; 
% sprintf('The time of finding the nearby voxel sequence is %s',tt)
% %% calculate the nearby voxel sequence
% 
% clear all;

end

