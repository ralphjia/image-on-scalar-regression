function [SRImgdata] = smoothResidual(mxR, Dimn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SVCM Project EXample1 June 30, 2011  LLK@CH %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input mxR n by N
% output smxR/SRImgdata n by N
% 
% aa = sprintf('%s/mxR.mat',resultPath);
% bb = sprintf('%s/smxR.mat',resultPath);

assert(ndims(mxR) == 2);
assert(size(mxR, 2) == prod(Dimn));
Kx = Dimn(1);
Ky = Dimn(2);
Kz = Dimn(3);

%% make Images Mask  
Mask = ones(Dimn);
id = find(Mask(:)~=0);
%% make Images Mask 

RImgdata = mxR;
SRImgdata = zeros(size(RImgdata));

%% Set up partition of the entire brain image
stride = 8;
Lx = 1:stride:Kx;
Ux = [Lx(2:end)-1, Kx];
Ly = 1:stride:Ky;
Uy = [Ly(2:end)-1, Ky];
Lz = 1:stride:Kz;
Uz = [Lz(2:end)-1, Kz];

Nx = Ux-Lx+1;
Ny = Uy-Ly+1;
Nz = Uz-Lz+1;
%% Set up partition of the entire brain image

NH = floor(max(10,max([Nx Ny Nz])/2));
HSeq = logspace(log10(1),log10(max([Nx Ny Nz])),NH); 

for ii=1:numel(Lx)
     
    XStartInd = Lx(ii);
    XEndInd = Ux(ii);
       
    for jj=1:numel(Ly)
        
        YStartInd = Ly(jj);
        YEndInd = Uy(jj);
             
        for kk=1:numel(Lz)

            % for each 8 x 8 x 8 consecutive cube:
            
            ZStartInd = Lz(kk);
            ZEndInd = Uz(kk);
            
            %% calculate voxel id with mask on in terms of the entire brain            
            % Not doing anything if masks == 1 for all voxels

            XYZCoord = zeros(Nx(ii)*Ny(jj)*Nz(kk),3);
            NonZeroCount = 0;
            for ZIndk = ZStartInd:ZEndInd
                for YIndj = YStartInd:YEndInd
                    for XIndi = XStartInd:XEndInd
                        if Mask(XIndi,YIndj,ZIndk)~=0
                           NonZeroCount = NonZeroCount + 1; 
                           XYZCoord(NonZeroCount,:) = [XIndi YIndj ZIndk];
                        end                        
                    end
                end
            end
            % 3D cube coordinates of the voxels in the cube
            XYZCoord = XYZCoord(1:NonZeroCount,:);
            % 1D cube coordinates of the voxels in the cube
            XYZInd = XYZCoord(:,1) + (XYZCoord(:,2)-1)*Dimn(1) + (XYZCoord(:,3)-1)*Dimn(1)*Dimn(2);
           
            % 1D imgage (i.e. global) coordinates of the voxels in the cube
            Tid = zeros(1,numel(XYZInd));
            mm=1;
            for ll=1:numel(XYZInd)   
                while id(mm)~= XYZInd(ll)
                    mm = mm+1;
                end
                Tid(ll) = mm;    
            end
            % Subset the cube out of the global imgage
            TRImgdata = RImgdata(:,Tid);
            
            
            %% Smooth individual residuals 

            if numel(XYZInd)>0
               [STRImgdata Opth,GCV] = siv(TRImgdata,XYZCoord,HSeq);% smooth residaul in mex code
               if Opth == 0
                   SRImgdata(:,Tid) = TRImgdata;
               else
                   SRImgdata(:,Tid) = STRImgdata;
               end
            end             
            
        end       
    end
end

end

