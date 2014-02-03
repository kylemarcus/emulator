function [T,P]=mytsearchn(x,tes,xi)
    sizex=size(x);
    sizexi=size(xi);
    sizetes=size(tes);
    if(~((numel(sizex)==2)&&(numel(sizexi)==2)&&(numel(sizetes)==2)&&...
         (sizex(2)==sizexi(2))&&(sizex(2)+1==sizetes(2))))
         error('bad data: wrong sizes');
    end
    Nxi=sizexi(1);
    Nx=sizex(1);
    Ndim=sizex(2);
    Ntes=sizetes(1);
    
    minx=min(min(x,[],1),min(xi,[],1));
    maxx=max(max(x,[],1),max(xi,[],1));
    invDX=((1+0.5/max(Nx,Nxi))*(maxx-minx)).^-1;

    x=x-repmat(minx,sizex(1),1);
    xi=xi-repmat(minx,sizexi(1),1);
    for idim=1:Ndim
        x(:,idim)=x(:,idim)*invDX(idim);
        xi(:,idim)=xi(:,idim)*invDX(idim);
    end

    Ndepth=ceil(log2(max(Ntes,Nxi))/Ndim)+1;
%     
%     [x,xBin,xIsort]=BinningSortId(x,Ndepth);
%     xIunsort(xIsort)=(1:Nx);
%     xIunsort=reshape(xIunsort,Nx,1);
%     tes=xIunsort(tes);
        
    xc=x';
    x=[x ones(Nx,1)]';
    xc=reshape(mean(reshape(xc(:,tes'),[Ndim Ndim+1 Ntes]),2),Ndim,Ntes)';    
    [xc,xcBin,xcIsort]=BinningSortId(xc,Ndepth);
    xc=[xc ones(Ntes,1)]';
    tes=tes(xcIsort,:);
    
    [xi,xiBin,xiIsort]=BinningSortId(xi,Ndepth);
    xiIunsort(xiIsort)=(1:Nxi);
    xiIunsort=reshape(xiIunsort,Nxi,1);
    xi=[xi ones(Nxi,1)]';
    
    MaxPtsNtes=ceil(Ntes*(Ndim+1)/Nx);
    Point2Simp=zeros(Nx,MaxPtsNtes);
    NsimpPerPoint=zeros(Nx,1);
    for ites=1:Ntes
        for inode=1:Ndim+1
            ipt=tes(ites,inode);
            NsimpPerPoint(ipt)=NsimpPerPoint(ipt)+1;
            Point2Simp(ipt,NsimpPerPoint(ipt))=ites;
        end
    end

    iafter=1;
    
    %might want to keep a list of the 10 (or so) last distinct simplices I 
    %found a point in and keep their J_L and J_U around so I can do a bunch
    %of good guesses cheaply before starting a search (i.e. having to
    %compute J_L and J_U for the new simplices I search)
    
    T=repmat(nan,Nxi,1);
    P=repmat(nan,Nxi,Ndim+1);
    iLastFoundSimp=0;
    NFoundPrevGuess=0;
    NFoundFirstGuess=0;
    NFoundSearch=0;
    NOutsideHull=0;
    for ixi=1:Nxi        
        iffound=0;
        %ifpointoutsideconvexhull=0;
        if(iLastFoundSimp)
            %xs(tes(iLastFoundSimp,:),:),
            p=CalcBaryCoord(xi(:,ixi),J_L_LastFound,J_U_LastFound);            
            if(all(p>=0))
                T(ixi)=xcIsort(iLastFoundSimp);
                P(ixi,:)=p';
                iffound=1;
                NFoundPrevGuess=NFoundPrevGuess+1;
            else
                J_L=J_L_LastFound;
                J_U=J_U_LastFound;
            end
        end

        if(~iffound)
            %find the simplex centers bounding (immediately before and 
            %after in the binning sorted list) the point we're looking for
            while((xiBin(ixi)>=xcBin(iafter))&&(iafter<Ntes))
                iafter=iafter+1;
            end            
            if(iafter>1)
                ibefore=iafter-1;
                %determine which of the bounding simplex centers is closer
                %to the current point
                d1=sum((xc(:,ibefore)-xi(:,ixi)).^2);
                d2=sum((xc(:,iafter )-xi(:,ixi)).^2);
                if(d1<=d2)
                    iGuessSimp=ibefore;
                else
                    iGuessSimp=iafter;
                end
            else
                iGuessSimp=iafter;
            end
                
            IGuessNode=tes(iGuessSimp,:);
            
            if(iGuessSimp==iLastFoundSimp)
                %already have p, J_L, J_U from above                                      
                NumRuleOutSimp=0;
                RuleOutSimp=[];
            else
                [p,J_L,J_U]=CalcBaryCoord(xi(:,ixi),x(:,IGuessNode));
                NumRuleOutSimp=1;
                RuleOutSimp=iLastFoundSimp;
            end
                
            if(all(p>=0))
                %hurray we found it with our first guess based on bounding
                %points
                iLastFoundSimp=iGuessSimp;
                J_L_LastFound=J_L;
                J_U_LastFound=J_U;
                T(ixi)=xcIsort(iLastFoundSimp);
                P(ixi,:)=p';
                NFoundFirstGuess=NFoundFirstGuess+1;
                %iffound=1;
            else
                %it wasn't in our first guess based on bounding simplex 
                %centers so we'll use the barycentric coordinates of this
                %first guess to choose the most likely neighbor, but just 
                %to make sure we don't check any simplex twice we'll keep
                %lists of simplices to rule out
                NumRuleOutSimp=NumRuleOutSimp+1;
                RuleOutSimp(NumRuleOutSimp)=iGuessSimp;
                p2=p;
                [p,isort]=sort(p,'descend');
                IGuessNode=IGuessNode(isort);
                nposp=find(p>=0,1,'last');                

                ifpointoutsideconvexhull=0;
                if((nposp==Ndim))
                    %Check to see if this simplex is on the boundary of the 
                    %convex hull and that point I'm looking for is on the
                    %otherside of that boundary.  If this is the (typical) 
                    %case this simplex will have Ndim positive nodes (a 
                    %face of the simplex) that aren't shared by any other 
                    %single simplex                   
                    NumFaceNeighList=NsimpPerPoint(IGuessNode(1))-1;
                    if(NumFaceNeighList)
                        FaceNeighList=setdiff(Point2Simp(IGuessNode(1),1:NsimpPerPoint(IGuessNode(1))),iGuessSimp);
                        FaceNeighListNodes=tes(FaceNeighList,:);
                        ifFaceNeigh=ones(numel(FaceNeighList),1);
                        for iposp=2:nposp
                            ifFaceNeigh=ifFaceNeigh&any(FaceNeighListNodes==IGuessNode(iposp),2);
                            NumFaceNeighList=sum(ifFaceNeigh);
                            if(NumFaceNeighList==0)
                                ifpointoutsideconvexhull=1;
                                NumSimpGuessList=0;
                                NOutsideHull=NOutsideHull+1;
                                %T(ixi)=-inf;
                                break;
                            end                        
                        end
                    end
                end
                    

                if(~ifpointoutsideconvexhull)                
                    %make a list of neighboring simplices to check (simplices
                    %that shared our last guessed simplex's most positive node
                    %(the node with the most positive barycentric coordinate) 
                    %and only consider the neigboring simplices that share all
                    %non-negative node and do not share exactly 1 negative node 
                    %rule out simplices we've already checked, and if there are
                    %multiple candidates choose the one whose center is closest
                    %to the current point
                    SimpGuessList=setdiff(Point2Simp(IGuessNode(1),1:NsimpPerPoint(IGuessNode(1))),RuleOutSimp);
                    NumSimpGuessList=numel(SimpGuessList);
                end
                while(NumSimpGuessList>0)
                    if(NumSimpGuessList==1)
                        iGuessSimp=SimpGuessList(1);
                    else
                        SimpGuessListNodes=tes(SimpGuessList,:);
                        SumPosp=repmat(p(1),NumSimpGuessList,1);
                        for iposp=2:nposp
                            SumPosp=SumPosp+any(SimpGuessListNodes==IGuessNode(iposp),2)*p(iposp);
                        end
                        [SumPosp,isort]=sort(SumPosp,'descend');
                        Nkeep=find(diff(SumPosp),1);
                        if(isempty(Nkeep))
                            Nkeep=NumSimpGuessList;
                        else
                            NumSimpGuessList=Nkeep;
                        end
                        i=isort(1:Nkeep);
                        SimpGuessList=SimpGuessList(i);
                        %might want to compare barycentric coordinates of
                        %SimpGuessList barycenters(in p3) to p2 here
                        SimpGuessListNodes=SimpGuessListNodes(i,:);
                        SumPosp=SumPosp(1);
                        if(SumPosp==p(1)&&(nposp>1))
                            %time to take a short cut
                            %p2=CalcBaryCoord(xi(:,ixi),J_L,J_U);
                            p3=CalcBaryCoord(xc(:,SimpGuessList),J_L,J_U);
                            d=(p3(1,:)-p2(1)).^2;
                            for idim=2:Ndim
                                d=d+(p3(idim,:)-p2(idim)).^2;
                            end
                            [garb,imin]=min(d);
                            SimpGuessList=SimpGuessList(imin);
                            NumSimpGuessList=1;
                        end
                        if(NumSimpGuessList==1)
                            iGuessSimp=SimpGuessList(1);
                        elseif(NumSimpGuessList>0)       
                            NumMissingNegp=zeros(NumSimpGuessList,1);
                            %save yadayada;
                            for inegp=nposp+1:Ndim+1
                                NumMissingNegp=NumMissingNegp+~any(SimpGuessListNodes==IGuessNode(inegp),2);
                            end
                            
                            [NumMissingNegp,isort]=sort(NumMissingNegp,'descend');
                            Nkeep=find(diff(NumMissingNegp),1);
                            if(isempty(Nkeep))
                                Nkeep=NumSimpGuessList;
                            else
                                NumSimpGuessList=Nkeep;
                            end                              
                            SimpGuessList=SimpGuessList(isort(1:Nkeep));
                            %NumMissingNegp=NumMissingNegp(1);
                        
                            if(Nkeep==1)
                                iGuessSimp=SimpGuessList(1);
                            elseif(NumSimpGuessList>1)
                                %save yadayada;
                                [garb,imin]=min(sum((xc(:,SimpGuessList)-repmat(xi(:,ixi),1,Nkeep)).^2));
                                iGuessSimp=SimpGuessList(imin);
                            end
                        end
                    end
                    if(NumSimpGuessList>0)
                        IGuessNode=tes(iGuessSimp,:);
                        [p,J_L,J_U]=CalcBaryCoord(xi(:,ixi),x(:,IGuessNode));
                        if(all(p>=0))
                            %hurray we found the simplex containing the
                            %current point
                            iLastFoundSimp=iGuessSimp;
                            J_L_LastFound=J_L;
                            J_U_LastFound=J_U;
                            T(ixi)=xcIsort(iLastFoundSimp);
                            P(ixi,:)=p';
                            NFoundSearch=NFoundSearch+1;
                            %iffound=1;
                            break;
                        end
                        NumRuleOutSimp=NumRuleOutSimp+1;
                        RuleOutSimp(NumRuleOutSimp)=iGuessSimp;
                        p2=p;
                        [p,isort]=sort(p,'descend');
                        IGuessNode=IGuessNode(isort);
                        nposp=find(p>=0,1,'last');                

                        ifpointoutsideconvexhull=0;
                        if((nposp==Ndim))
                            %Check to see if this simplex is on the boundary of the 
                            %convex hull and that point I'm looking for is on the
                            %otherside of that boundary.  If this is the (typical) 
                            %case this simplex will have Ndim positive nodes (a 
                            %face of the simplex) that aren't shared by any other 
                            %single simplex                   
                            NumFaceNeighList=NsimpPerPoint(IGuessNode(1))-1;
                            if(NumFaceNeighList)
                                FaceNeighList=setdiff(Point2Simp(IGuessNode(1),1:NsimpPerPoint(IGuessNode(1))),iGuessSimp);
                                FaceNeighListNodes=tes(FaceNeighList,:);
                                ifFaceNeigh=ones(numel(FaceNeighList),1);
                                for iposp=2:nposp
                                    ifFaceNeigh=ifFaceNeigh&any(FaceNeighListNodes==IGuessNode(iposp),2);
                                    NumFaceNeighList=sum(ifFaceNeigh);
                                    if(NumFaceNeighList==0)
                                        ifpointoutsideconvexhull=1;
                                        NumSimpGuessList=0;
                                        NOutsideHull=NOutsideHull+1;
                                        %T(ixi)=-inf;
                                        break;
                                    end                        
                                end
                            end
                        end
                    
                        if(~ifpointoutsideconvexhull)
                            %make a list of neighboring simplices to check (simplices
                            %that shared our last guessed simplex's most positive node
                            %(the node with the most positive barycentric coordinate) 
                            %and only consider the neigboring simplices that share all
                            %non-negative node and do not share exactly 1 negative node 
                            %rule out simplices we've already checked, and if there are
                            %multiple candidates choose the one whose center is closest
                            %to the current point
                            SimpGuessList=setdiff(Point2Simp(IGuessNode(1),1:NsimpPerPoint(IGuessNode(1))),RuleOutSimp);
                            NumSimpGuessList=numel(SimpGuessList);
                        end        
                    end    
                end
                

 
            end
        end
    end     
          
    %NFound=NFoundPrevGuess+NFoundFirstGuess+NFoundSearch
    %NFoundPrevGuess
    %NFoundFirstGuess
    %NFoundSearch
    %NOutsideHull
    
    T=T(xiIunsort);
    P=P(xiIunsort,:);
    %save DEBUGME; 
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Coord,Bins,Isort]=BinningSortId(Coord,Ndepth,varargin)

[Npts,Ndim]=size(Coord);
if(nargin<3)
    Isort=(1:Npts)';
else
    Isort=reshape(varargin{1},Npts,1);
end
    
Coord=2*Coord;
temp=floor(Coord);

[Bins,isort]=sort(temp*2.^(0:Ndim-1)'); % sum(temp,2)],[-1]);
Isort=Isort(isort);
Bins=Bins(:,1);
temp=temp(isort,:);
Coord=Coord(isort,:);
ScaleBins=0.5^Ndim;

if((Ndepth<=1)&&all(Coord(1,:)==Coord(Npts,:)))
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
    [coord,bins,isort]=BinningSortId(Coord(i,:),Ndepthm1,Isort(i));
    Coord(i,:)=0.5*(coord+temp(i,:));
    Bins(i)=ScaleBins*(Bins(i)+bins);
    Isort(i)=isort;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,J_L,J_U]=CalcBaryCoord(x,varargin)

if(nargin==2)
    [J_L,J_U]=lu(varargin{1});
else
    J_L=varargin{1};
    J_U=varargin{2};
end
%J_L is either a lower triangular or permuted lower a.k.a "psychologically
%lower" triangular matrix so can't use linsolve with it
opts.UT = true;
P=linsolve(J_U,(J_L\x),opts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ifOutside,p,IGuessNode]=CheckIfOutsideConvexHull(iGuessSimp,p,tes)
% 
% 
% 
%                         ifpointoutsideconvexhull=0;
%                         if((nposp==Ndim))
%                             %Check to see if this simplex is on the boundary of the 
%                             %convex hull and that point I'm looking for is on the
%                             %otherside of that boundary.  If this is the (typical) 
%                             %case this simplex will have Ndim positive nodes (a 
%                             %face of the simplex) that aren't shared by any other 
%                             %single simplex                   
%                             NumFaceNeighList=NsimpPerPoint(IGuessNode(1))-1;
%                             if(NumFaceNeighList)
%                                 FaceNeighList=setdiff(Point2Simp(IGuessNode(1),1:NsimpPerPoint(IGuessNode(1))),iGuessSimp);
%                                 FaceNeighListNodes=tes(FaceNeighList,:);
%                                 ifFaceNeigh=ones(numel(FaceNeighList),1);
%                                 for iposp=2:nposp
%                                     ifFaceNeigh=ifFaceNeigh&any(FaceNeighListNodes==IGuessNode(iposp),2);
%                                     NumFaceNeighList=sum(ifFaceNeigh);
%                                     if(NumFaceNeighList==0)
%                                         ifpointoutsideconvexhull=1;
%                                         NumSimpGuessList=0;
%                                         NOutsideHull=NOutsideHull+1;
%                                         break;
%                                     end                        
%                                 end
%                             end
%                         end
