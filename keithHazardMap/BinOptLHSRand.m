function coord=BinOptLHSRand(NdimWant,NptsWant)
%this function will generate binning optimal symmetric latin hypercube
%designs for Ndim=1,2,4, or 8 dimensions, with a number of points that is
%equal to 2*Ndim*2^p  where p is a non negative integer.  The math requires
%a number of operations of order O(Npts*log(Npts)) however because MATLAB's
%built in quicksort is a stable sort (relative order of entries with the
%same value does not change) the cost is O(Npts/Ndim*(log(Npts))^2) instead
%the function returns sample points whose coordinates lie between 0 and 1.

if(ischar(NdimWant)||ischar(NptsWant))
   error('Do not call this function from outside MATLAB'); 
end
rand('twister',5489);

    [Ndim,Nref,Nbits]=PickNdimNrefNbits(NdimWant,NptsWant);
    %NptsTotal=2*Ndim*2^Nref;
    coord=GenOrthoBaseLhsByTiling(Ndim);    
    [Npts,Ndim]=size(coord);
    signs=sign(coord(1:Ndim,:));

    %external version
    %[orient,priority]=GenSortedDistinctOrientations(Ndim);
    %Orient=2*(dec2bin(orient,Ndim)=='1')-1;
    %Orient=Orient(:,Ndim:-1:1)

    %internal version
    [Orient,priority]=GenSortedDistinctOrientations(Ndim);

    coord=0.5*(coord/Npts+1);
    for iref=1:Nref
        ndepth=ceil(log2(Npts)/Ndim);
        [coord,bin]=BinningSortId(coord,ndepth);
        coord=AssignRandBinOptOct(coord,Orient,signs);    
        coord=0.5*([coord; -coord]+1);
        Npts=Npts*2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ndim,Nref,Nbits]=PickNdimNrefNbits(Ndimwant,Nptswant)
%This function corrects bad inputs provided by the user.
%Bad means any of the following
%*the number of dimensions is not a power of 2
%*the number of is not (2*Ndim)*2^Nref (integer Nref)
%*the number of desired points is > 2^64 (more than 64 bits would be
% required to uniquely represent the position of any point).
%
%This function also determines the number of bits required to represent
%location (choosing between 8, 16, 32, 64) and the number of refinements
%(doubling the current number of points) needed to generate the total
%number of points.

    Ndim=2^ceil(log2(Ndimwant));
    Nrefwant=ceil(log2(Nptswant/(2*Ndim)));
    Nptswant2=2*Ndim*2^Nrefwant;
    Nbits=min(64,8*2^ceil(max(0,log2(Nptswant2)/8-1)));
    Nref=min(Nrefwant,Nbits-(1+log2(Ndim)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coord=GenOrthoBaseLhsByTiling(Ndim)
%this function generates an orthogonal LHS design of Ndim points in
%Ndim dimensions and then "balances" the design by adding reflections 
%through zero of original points.  It generates the orthogonal design
%by "tiling" separate patterns for the coordinates' magnitudes and 
%signs and combining the 2 patterns in a 1 for 1 fashion.  In order for
%the pattern to work Ndim must be a power of 2, so if Ndim is not a 
%power of 2 this function increases Ndim to the next larger power of 2.
%The output of this function is a Npts by Ndim array of coordinates where 
%(Npts=2*Ndim and) each row of coord is a point in Ndim dimensional space.
%half, all, or none of the signs on each point are positive, and every
%quadrant of every 2 dimensional pairing of dimensions contains the same
%number of points.  Because this function generates the orthogonal 
%coordinates by tiling simple patterns (rather than complicated
%calculations), it is very fast.

    Ntimes2tile=ceil(log2(Ndim));
    Ndim2=2^Ntimes2tile;

    %tile the magnitudes of the coordinates and signs separately then
    %combine
    coord=zeros(Ndim2);
    maxnum=2*Ndim2-1;
    coord(1,:)=1:2:maxnum;
    signs=1;
    for Itile=1:Ntimes2tile
        tilesize=2^Itile;
        i1=1:tilesize/2;
        i2=tilesize/2+1:tilesize;
        jj1=1:tilesize;
        jj2=tilesize:-1:1;
        for jtile=1:Ndim2/tilesize
            coord(i2,jj2+(jtile-1)*tilesize)=coord(i1,jj1+(jtile-1)*tilesize);
        end
        signs=[signs -signs; signs signs]; %2^0.5* rotation (by pi/4)
        %matrix sines and cosines are orthogonal functions, this is a
        %Hadamard matrix
    end

    coord=signs.*coord;
    coord=[coord; -coord];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Orient,priority]=GenSortedDistinctOrientations(Ndim)
%generate a list of non-overlaping orientations arranged in a maximally 
%spaced order, this is an analytical solution for Ndim=1,2,4,8
    if((Ndim==1)||(Ndim==2))
        Orient=ones(1,Ndim); priority=0;
    elseif(Ndim==4)
        Orient=ones(2,Ndim); Orient(2,1)=-1;
        Orient(2,:)=-Orient(2,:); priority=1;
    elseif(Ndim==8)
        Orient=repmat(ones(Ndim)-2*eye(Ndim),2,1);
        Orient(1:8,1)=-Orient(1:8,1);
        Orient(2:9,:)=-Orient(2:9,:);
        priority=reshape(repmat([2 1],8,1),16,1);        
    else
        error('GenSortedDistinctOrientations: Currently only works for Ndim=2^N, N=0:3');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Coord,Bins]=BinningSortId(Coord,Ndepth)
%Ndepth is the depth of the binning (recursive division of hypercube edges)
%a binning quicksort that computes bin id's in the process

    [Npts,Ndim]=size(Coord);

    Coord=2*Coord;
    temp=floor(Coord);

    [Bins,isort]=sort(temp*2.^(0:Ndim-1)'); % sum(temp,2)],[-1]);

    Bins=Bins(:,1);
    temp=temp(isort,:);
    Coord=Coord(isort,:);
    ScaleBins=0.5^Ndim;

    if(Ndepth<=1)
        Coord=0.5*Coord;
        Bins=ScaleBins*Bins;
        return;
    end
    Ndepthm1=Ndepth-1;
    Coord=Coord-temp;

    i=find(diff(Bins));
    NBinNotEmpty=numel(i)+1;
    iStartEnd=[[1; i+1] [i; Npts]];

    for iBinNotEmpty=1:NBinNotEmpty
        i=iStartEnd(iBinNotEmpty,1):iStartEnd(iBinNotEmpty,2);
        [coord,bins]=BinningSortId(Coord(i,:),Ndepthm1);
        Coord(i,:)=0.5*(coord+temp(i,:));
        Bins(i)=ScaleBins*(Bins(i)+bins);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coord=AssignRandBinOptOct(coord,Orient,signs)

    [Npts,Ndim]=size(coord);
    if(Ndim==1)
        return;
    elseif(Ndim==2)
        NMaxSpaceOrient=1;
    elseif(Ndim==4)
        NMaxSpaceOrient=2;
    else
        NMaxSpaceOrient=8;
    end
    NOrient=size(Orient,1);

    if(Npts<2^Ndim)
        IOrtho=GenRandIOrtho(Ndim);
        IOrient=GenRandIOrient(Ndim,NOrient,NMaxSpaceOrient);
        i=1:Npts;
        coord=coord.*signs(IOrtho(i),:).*Orient(IOrient(i),:);
    else
        TwoToNdim=2^Ndim;
        i=1:TwoToNdim;
        Nbatch=Npts/TwoToNdim;
        for ibatch=1:Nbatch
            IOrtho=GenRandIOrtho(Ndim);
            IOrient=GenRandIOrient(Ndim,NOrient,NMaxSpaceOrient);
            %save DEBUGME;
            coord(i,:)=coord(i,:).*signs(IOrtho,:).*Orient(IOrient,:);
            i=i+TwoToNdim;
        end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IOrtho=GenRandIOrtho(Ndim)
    HalfTwoToNdim=2^(Ndim-1);
    Nhalf=HalfTwoToNdim/Ndim;
    IOrtho=zeros(Ndim,2*Nhalf);
    for i=1:Nhalf
        IOrtho(:,i)=randperm(Ndim)';
    end
    i=1:HalfTwoToNdim;
    IOrtho(2*HalfTwoToNdim+1-i)=IOrtho(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IOrient=GenRandIOrient(Ndim,NOrient,NMaxSpaceOrient)
    IOrient=1:NOrient;
    yada=NOrient/NMaxSpaceOrient-1;
    for iyada=0:yada
        i1=iyada*NMaxSpaceOrient+(1:NMaxSpaceOrient);
        i2=iyada*NMaxSpaceOrient+randperm(NMaxSpaceOrient);
        IOrient(i2)=IOrient(i1);
    end
    IOrient=repmat([IOrient IOrient(NOrient:-1:1)],Ndim,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%