function build_mini_emulator(samplenumber)
    rand('twister',5489);
    iftoc=1;
    if(iftoc)
        tic;
    end
    ifplot=0; %no plots currently possible, 
    ifwrite=1;
    if(ischar(samplenumber))
        samplenumber=str2double(samplenumber);
    end
    
    fid=fopen(sprintf('build_mini_pwem_meta.%06g',samplenumber),'r');
    Nskip=sscanf(fgets(fid),'additional file format lines=%g',1);
    for i=1:Nskip
        fgets(fid);
    end

    Ndiminmacro=sscanf(fgets(fid),'%g',1);
    lminmacro=sscanf(fgets(fid),'%g',[1 Ndiminmacro]);
    Nneighmacro=sscanf(fgets(fid),'%g',1);
    iyneighmacro=fscanf(fid,'(%g) %g\n',[2 Nneighmacro])';
    fclose(fid);

    %partly hardcoded for Ndimout=1;
    Ndimout=1;
    Nspacedim=2;
    Ndimin=Ndiminmacro+Nspacedim;
    Idimin=1:Ndimin;
    Idimout=Ndimin+1;
    idimout=1;
    
    ysy=[];
    YMACRO=[];
    ieastnorthh=[1 2 Idimout];
    Nhop=2;
    for ineighmacro=1:Nneighmacro
        jyneighmacro=iyneighmacro(ineighmacro,2);
        filename=sprintf('down_sampled_data.%06g',jyneighmacro);
        fid=fopen(filename,'r');
        
        Nskip=sscanf(fgets(fid),'additional file format lines=%g',1);
        for i=1:Nskip
            fgets(fid);
        end
        
        %fprintf(fid,'%%Ndiminmacro=4: of macro-emulator input dimensions (log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
        %fprintf(fid,'%%ymacro=(log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg]): values of "uncertain inputs" for the particular simulation: [Ndiminmacro] doubles\n');
        %fprintf(fid,'%%Ny Nycover: #of point and #of points needed to cover the flow: [2] integers\n');
        %fprintf(fid,'%%[Nycover] lines of point info: each contains {{(iy): pt id: [1] integer},{y=(east, north, h) [3] doubles},{Nneigh: # of pts in neighborhood: [1] integer},{iyneigh: indices of pts in neighborhood: [Nneigh] integers},{neighhops: #number of hops each pt in neighborhood is away from iy: [Nneigh] integers}}\n');
        %fprintf(fid,'%%[Ny-Nycover] lines of limited point info: each contains {{(iy): pt id: [1] integer},{y=(east,north): we know h=0 for these: [2] doubles}}\n');

        %fprintf(fid,'%%{Ntri Ntricover}: #of triangles in east-north tessellation and the #of triangles needed to cover the flow: [2] integers\n');
        %fprintf(fid,'%%tri: triangle node indices (into y): [Ntri 3] array of doubles\n');

        Ndiminmacrotmp=sscanf(fgets(fid),'%g',1);
        if(Ndiminmacrotmp~=Ndiminmacro)
            error('bldminibaddim','build_mini_emulator(): iyneighmacro(%g)=%g: Ndimimacro(iyneighmacro)=%g ~= Ndimimacro=%g',...
                ineighmacro,jyneighmacro,Ndiminmacrotmp,Ndiminmacro);
        end
        ymacrotmp=sscanf(fgets(fid),'%g',[1 Ndiminmacro]);
        YMACRO=[YMACRO; ymacrotmp];
        
        if(ineighmacro==1)
            yada=sscanf(fgets(fid),'%g %g',2);
            Ny0T=yada(1);
            Ny0N=yada(2);
            Ny0E=Ny0T-Ny0N;

            Ny=Ny0T; 
            ysy=[zeros(Ny0T,2) repmat(ymacrotmp,Ny0T,1) zeros(Ny0T,1)];
            
            pwem=repmat(struct(...
                'Nneigh0',0,'iy0neigh',[],'i0hop',[],...
                'Nneigh' ,0,'iyneigh' ,[],'ihop' ,[],...
                'lmin',[],'lmax',[],'theta',[],'beta',[],'sigma',inf,...
                'crit',inf,'Rinveps',[]),Ny0N,1);
            
            for iy0=1:Ny0N
                current_input_line=fgets(fid);
                yada=sscanf(current_input_line,'(%g) %g %g %g %g',5);
                Nneigh0=yada(5);            
                
                yada=sscanf(current_input_line,sprintf('(%%g) %%g %%g %%g %%g%s',repmat(' %g',1,2*Nneigh0)),5+2*Nneigh0);
                ysy(iy0,ieastnorthh)=yada(2:4);
                
                pwem(iy0).Nneigh0=Nneigh0;
                pwem(iy0).iy0neigh=yada(5+(1:Nneigh0));
                pwem(iy0).i0hop=yada(5+Nneigh0+(1:Nneigh0));
                
                pwem(iy0).Nneigh=pwem(iy0).Nneigh0;
                pwem(iy0).iyneigh=pwem(iy0).iy0neigh;
                pwem(iy0).ihop=[ones(Nneigh0,1) zeros(Nneigh0,1) pwem(iy0).i0hop];
            end
            
            yada=fscanf(fid,'(%g) %g %g\n',[3 Ny0E])';
            ysy(Ny0N+1:Ny0T,1:2)=yada(:,2:3);
            
            for iy0=1:Ny0N
                %set the minimum and maximum correlation lengths to guess 
                %between. This limits the rate at which adaptive sampling
                %would be able to _tightly_place_ points                
                iy0hop1=pwem(iy0).iy0neigh(find(pwem(iy0).i0hop<=1));              
                if(1)
                    %"new and improved" logic, in the process of being
                    %tested...
                    %lmin=1/4 the radius of the circle with the same 
                    %area as the sum of the 1 hop triangles (only 2 new
                    %dimensions are east and north and they have the same
                    %scale). The radius of the 1 hop circle equals the 
                    %average distance to adjacent points/simulations.  We
                    %want a "completely unsure" micro-emulator to say it 
                    %is "most unsure" about the point 1/2 this distance 
                    %away so say 1/2 radius = 2 correlation lengths 
                    %(standard deviations) 
                    %=> minium correlation length, lmin, = radius/4
                    iy0hop1convhull=iy0hop1(convhull(ysy(iy0hop1,1),ysy(iy0hop1,2)));
                    pwem(iy0).lmin=[[0.25 0.25]*sqrt(polyarea(...
                        ysy(iy0hop1convhull,1),ysy(iy0hop1convhull,2)...
                        )/pi) lminmacro];
                else
                    %For each dimension use 1/4th the average distance to 
                    %1 hop neighbors as the minimum correlation length this 
                    %micro-emulator.  This means that a "locally completely
                    %unsure ensemble emulator" will say it is most 
                    %uncertain at a point approximately at the center of
                    %the simplex.  
                    pwem(iy0).lmin=[0.25*sqrt(mean(...
                        (ysy(iy0hop1,1:Nspacedim)-...
                         repmat(ysy(iy0,1:Nspacedim),[length(iy0hop1) 1])...
                        ).^2)*Nspacedim) lminmacro];
                end
                %we want the maximum correlation length to be 2 hops of 
                %distance, lmin=(1 hop dist)/4 => lmax=8*lmin, data points
                %further away than 4*(1 hop dist) from this point (iy0) 
                %_MAY_ ALSO be trustworthy, but with a gaussian error model
                %that would mean that the same information would be 
                %contained in the closer data, so we don't need to go 
                %farther than that.
                pwem(iy0).lmax=8*pwem(iy0).lmin;                               
            end

            %read in the spatial (east,north) tessellation for the primary
            %simulation
            yada=sscanf(fgets(fid),'%g %g',2);
            Ntess0T=yada(1);
            Ntess0N=yada(2);
            tess0=fscanf(fid,'%g',[3 Ntess0T])';
            
            Iy0N=1:Ny0N;            
        else
            %read in data points from macro neighbor simulations
            
            Nyb4=Ny; %total number of data points that have come before now
   
            yada=sscanf(fgets(fid),'%g %g',2); %the number of new data points
            %being added from this simulation.
            NynewT=yada(1);
            NynewN=yada(2);
            NynewE=NynewT-NynewN;

            
            ysynew=[zeros(NynewT,2) repmat(ymacrotmp,NynewT,1) zeros(NynewT,1)];
            %the new data points being added from this simulation

            %save DEBUMEYADA;
            Nhopspacekeep=Nhop-iyneighmacro(ineighmacro,1); %define the full 
            %dimensional neighborhood of each primary simulation data point
            %to be all data points whose total (add) number of hops in the
            %spatial (east,north) PLUS (add) macro dimension from the 
            %primary data point at the center of the micro-emulator is less
            %than or equal to Nhop hops, therefore the number of allowed
            %spatial hops in this non-primary simulation is Nhopspacekeep
            %I did this so there would be fewer points in each micro
            %emulator, to speed up the building the micro emulator, it was
            %taking far too long before I did this, in one case it took 7
            %hours for a single mini-emulator with about 1800
            %micro-emulators
            
            pwemnew=repmat(...
                struct('Nneigh',0,'iyneigh',[],'ihop',[]),...
                NynewN,1);
            
            for iynew=1:NynewN
                current_input_line=fgets(fid);
                yada=sscanf(current_input_line,'(%g) %g %g %g %g',5);
                Nneigh=yada(5);            
                yada=sscanf(current_input_line,sprintf('(%%g) %%g %%g %%g %%g%s',repmat(' %g',1,2*Nneigh)),5+2*Nneigh);
                ysynew(iynew,ieastnorthh)=yada(2:4);                
                iyneigh=Nyb4+yada(5+(1:Nneigh)); %added Nyb4 to make it the indice into entire list of ysy
                ihop=yada(5+Nneigh+(1:Nneigh));
                ikeep=find(ihop<=Nhopspacekeep); %theses are the kept spatial hops 
                pwemnew(iynew).Nneigh=numel(ikeep);
                pwemnew(iynew).iyneigh=iyneigh(ikeep);
                pwemnew(iynew).ihop=ihop(ikeep);
            end
            yada=fscanf(fid,'(%g) %g %g\n',[3 NynewE])';
            ysynew(NynewN+1:NynewT,1:2)=yada(:,2:3);            
            ysy=[ysy; ysynew];
            Ny=Nyb4+NynewT;
            
            yada=sscanf(fgets(fid),'%g %g',2);
            NtessnewT=yada(1);
            NtessnewN=yada(2);
            tessnew=fscanf(fid,'%g',[3 NtessnewT])';
            
            Iynearestnew=dsearch(ysynew(:,1),ysynew(:,2),tessnew,ysy(Iy0N,1),ysy(Iy0N,2));
            %dsearchn(ysynew(:,1:2),tessnew,ysy(Iy0N,1:2));
            IynewN=1:NynewN;
            for iy0=1:Ny0N
               iynew=Iynearestnew(iy0);
               if(iynew>NynewN)
                  [garb,iynew]=min(sum((repmat(ysy(iy0,1:2),NynewN,1)-ysynew(IynewN,1:2)).^2,2));                   
               end
               pwem(iy0).Nneigh=pwem(iy0).Nneigh+pwemnew(iynew).Nneigh;
               pwem(iy0).iyneigh=[pwem(iy0).iyneigh; pwemnew(iynew).iyneigh];
               pwem(iy0).ihop=[pwem(iy0).ihop;...
                   repmat([ineighmacro iyneighmacro(ineighmacro,1)],pwemnew(iynew).Nneigh,1) pwemnew(iynew).ihop];
            end
        end
        fclose(fid);
    end %for ineighmacro=1:Nneighmacro
    
    if(Ny~=size(ysy,1))
        save DEBUG_BUILD_MINI_EMULATOR0;
        disp(sprintf('ERROR: Ny=%g ~= size(ysy,1)=%g: saved workspace to DEBUG_BUILD_MINI_EMULATOR0.mat\n',Ny,size(ysy,1)));
        bob; %cause the code to crash
    end

    clear jyneighmacro filename Nskip ymacrotmp pwemnew ysynew iynew Nynew...
        Iynearestnew tessnew Ntessnew ihop ikeep iyneigh Nhopspacekeep...
        Nyb4 yada;
    
    %save CLEANUP1;
    

    
    if(iftoc)
        disp(sprintf('done loading downsampled data... Ny0N=%g Ny0T=%g Ny=%g at t=%g [sec]',Ny0N,Ny0T,Ny,toc));
    end
          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %build and write the micro emulators
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ifwrite)
        Ndiminspatial=Ndimin-Ndiminmacro;
        
        fid=fopen(sprintf('mini_emulator.%06d',samplenumber),'w');
        fprintf(fid,'additional file format lines=%g\n',21);
        
        %this data can be printed here
        fprintf(fid,'%%iymacro: index/id of the simulation at the center of this mini-emulator, this should be the same as the filename suffix: [1] integer\n');
        fprintf(fid,'%%Ndiminspatial=2: #number of spatial input dimension (east,north): [1] integer\n');
        fprintf(fid,'%%Ndimin=6: # of input dimensions (east,north,log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
        fprintf(fid,'%%Ny0T: #of data points (y0) in this mini-emulator''s central simulation: [1] integer\n');
        fprintf(fid,'%%Ny0N: #of micro-emulators in this mini-emulator: [1] integer\n');
        fprintf(fid,'%%Ny: #of data points in this mini-emulator: [1] integer\n');        
        fprintf(fid,'%%{{y=(east,north,log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg])},{sy=h}} [Ndimin+1] doubles, there are [Ny] of these lines, the first [Ny0N] of them are the centers of the micro-emulators\n');
        fprintf(fid,'%%NsimpT: # of simplices in y0 tesselation: [1] integer\n');
        fprintf(fid,'%%NsimpN: # of simplices in y0 tesselation that have at least 1 micro-emulator as a node: [1] integer\n');
        fprintf(fid,'%%tess0: [NsimpT Ndimin+1] array of integers/data-point-indices\n');
        
        %this data must be printed in the following loop
        fprintf(fid,'%%then there will be [Ny0N] micro-emulator "structures," each of which contains the following...\n');
        fprintf(fid,'%%(iy0): index of micro emulator, should also be the first entry in iy0neigh and iyneigh: [1] integer\n');
        fprintf(fid,'%%lmax: maximum plausible correlation lengths: [Ndimin] doubles\n');
        fprintf(fid,'%%lmin: minimum plausible correlation lengths: [Ndimin] doubles\n');
        fprintf(fid,'%%theta: roughness parameters: [Ndimin] doubles\n');
        fprintf(fid,'%%sigma2: unadjusted variance: [1] double\n');
        fprintf(fid,'%%beta: least squares coefficients: [Ndimin+1] doubles\n');
        fprintf(fid,'%%Nneigh0: #of y0 data points in the neighborhood of this micro-emulator: [1] integer\n');
        fprintf(fid,'%%{{iy0neigh: index of y0 neighbor: [1] integer},{i0hop: #of y0 hops from center of micro-emulator: [1] integer}}\n');
        fprintf(fid,'%%Nneigh: #of y data points in the neighborhood of this micro-emulator: [1] integer\n');
        fprintf(fid,'%%{{iyneigh: index of y neighbor: [1] integer},{ihop: #of (y0,ymacro) hops from center of micro-emulator: [3] integers},{Rinveps: R^-1*eps (eps = sy minus unadjusted micro-emulator mean): [1] double}}\n');
        
        fprintf(fid,'%g\n',samplenumber,Ndiminspatial,Ndimin,Ny0T,Ny0N,Ny);
        yada=sprintf('%%.2f %%.2f%s\n',repmat(' %.10g',1,Ndiminmacro+1));
        fprintf(fid,yada,ysy');
        fprintf(fid,'%g\n',Ntess0T,Ntess0N);
        fprintf(fid,'%g %g %g\n',tess0');        
        
        thetaformat=sprintf('%%g%s\n',repmat(' %g',1,Ndimin-1));
        betaformat=sprintf('%%g%s\n',repmat(' %g',1,Ndimin));
    end
    
    idiminzeroh=3:Ndimin; 
    for iy0=1:Ny0N
    %for iy0=399:401;
        if(iftoc)
            disp(sprintf('iy0=%g/%g building micro-emulator at t=%g [sec]',iy0,Ny0N,toc));
        end
        pwemthis=pwem(iy0);
        
        Nyn=pwemthis.Nneigh;
        ysyn=ysy(pwemthis.iyneigh,:);

        ONES=ones(Nyn,1);

        if((ysyn(1,Idimout)==0))
            %this "hack" means that for any data point with zero pile
            %height the coefficients for the height and the east and north
            %direction slopes will be zero
            A=ysyn(:,idiminzeroh)-repmat(ysyn(1,idiminzeroh),[Nyn 1]);
            pwemthis.beta=[0; 0; 0; (A'*A)\(A'*ysyn(:,Idimout))];
            epsn=ysyn(:,Idimout);
        else
            A=[ONES ysyn(:,Idimin)-repmat(ysyn(1,Idimin),[Nyn 1])];
            pwemthis.beta=(A'*A)\(A'*ysyn(:,Idimout));
            epsn=ysyn(:,Idimout)-max(A*pwemthis.beta,0);
        end    

        Zn=zeros(Nyn^2,Ndimin);
        for idimin=1:Ndimin
            x=repmat(ysyn(:,idimin),[1,Nyn]);
            Zn(:,idimin)=-reshape((x-x').^2,[Nyn^2 1]);
        end
        Nfree=max(Nyn-(Ndimin+1),1);

        if(0)
            %this does a non-trivial bit better with the adjusted mean but
            %much worse at the adjusted variance 
            alpha=1/Nyn;
            b=(alpha^2 + (1-alpha)^2-Nyn/Nfree)/(alpha*(1-alpha));
            minrcondR=1/((-b+(b.^2-4).^0.5)/2);
        else
            %minrcondR=singularvalue();
            minrcondR=1/Nfree;
        end
    
        CorLenBnds=struct('Ndimin',Ndimin,'lmax',pwemthis.lmax,...
            'lmin',pwemthis.lmin,'Pmax',log2(pwemthis.lmax./pwemthis.lmin));
        
        NThetaGuesses=5;
                
        for itheta=1:NThetaGuesses

            [crit,sigma,theta,r]=GuessGoodThetaSigma(4,Nfree,epsn,...
                Zn,minrcondR,CorLenBnds,ysyn);
            %itheta
            %crit 
            %sigma
            %theta
            if((sigma>=0)&&(crit<=pwemthis.crit))
                pwemthis.crit =crit;
                R             =r;
                pwemthis.sigma=sigma;        
                pwemthis.theta=theta;
            end
        end
        %save DEBUGMEYADA2;
        pwemthis.Rinveps=R\epsn;
    
        if(ifwrite)
            fprintf(fid,'(%g)\n',iy0);
            fprintf(fid,thetaformat,pwemthis.lmax);
            fprintf(fid,thetaformat,pwemthis.lmin);
            fprintf(fid,thetaformat,pwemthis.theta);
            fprintf(fid,'%.10g\n',pwemthis.sigma);
            fprintf(fid,betaformat,pwemthis.beta);
            fprintf(fid,'%g\n',pwemthis.Nneigh0);
            fprintf(fid,'%g %g\n',[pwemthis.iy0neigh pwemthis.i0hop]');
            fprintf(fid,'%g\n',pwemthis.Nneigh);            
            fprintf(fid,'%g %g %g %g %g\n',[pwemthis.iyneigh pwemthis.ihop pwemthis.Rinveps]');
        end
            
        %pwemthis.Vary=diag(emulator.sigma);
        %pwemthis.invVary=diag(emulator.sigma.^-1);
        %pwemthis.gamma=eye(Ndimout);
    end
    %bob;

    if(ifwrite)
        fclose(fid);
    end

    %emulator=struct('Ndimin',Ndimin,'Ndimout',Ndimout,...
    %    'Ny',Ny,'y',ysy(:,Idimin),'sy',ysy(:,Idimout),'fy',fy,'tess0',tess0,'pwem',pwem);
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function singularval=singularvalue()
    singularval=eps;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fy=falsey(fy)
  %"move" points nearest corners of the hypercube to corners of the
  %hypercube
  Ny=size(fy,1);
  Ndimin=size(fy,2);
  
  Nnodes=2^Ndimin;
  %nodes=zeros(max(Nnodes+1,Nmaxnodes),Ndimin); 
  inode=zeros(1,Ndimin); %hyper indices of hypercube corners
  inode(1)=-1;
  for i=1:Nnodes
    inode(1)=inode(1)+1;
    
    for idimin=1:Ndimin-1
      if(inode(idimin)>1)
        inode(idimin)=0;
        inode(idimin+1)=inode(idimin+1)+1;
      else
        break
      end
    end
      
    %nodes(i,:)=2*inode-1; 
    node=2*inode-1; 
    [garb,iy]=min(sum((fy-repmat(node,[Ny 1])).^2,2));
    fy(iy,:)=node;
  end
  

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ALL code below this line is used to produce the error model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [CRIT,SIGMA,THETA,R]=...
    GuessGoodThetaSigma(Nnewt,Nfree,eps,...
        Z,minrcondR,CorLenBnds,ysy)
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:illConditionedMatrix');

    SmallNearlySingular=singularvalue();

    Nsample=size(eps,1);
    %Ndimin=size(Z,2); %don't need this in this function
    
    MinTheta=0.5.*CorLenBnds.lmax.^-2;
    %MaxTheta=0.5.*CorLenBnds.lmin.^-2; 
        
    THETA=rand_theta(CorLenBnds); %guess a random theta
    %save DEBUGMEYADA1;
    [THETA,R,rcondR,eigvect,eigval,d1rcondRdti,Q1,Q2]=...
        increase_rcondR_to_minrcondR(Z,minrcondR,THETA);
        
    %save DEBUGMEYADA2;
    %I=eye(Nsample);
    %[Rinv,ifsing]=chol(R);
    %if(ifsing)
    %    save DEBUG_BUILD_MINI_EMULATOR1;
    %    disp('ERROR: singular R!!!! Workspace has been saved to DEBUG_BUILD_MINI_EMULATOR1.mat\n');
    %    bob; %cause the code to crash
    %end
    %
    %Rinv=Rinv\(Rinv'\I);

    if(~any(eps))
        SIGMA=0;
        CRIT=0;
        return;
    end
    
    if(any(THETA<MinTheta))
    %if(any(THETA<MinTheta)||any(THETA>MaxTheta))
    %if(any(THETA<MinTheta)||any(THETA>MaxTheta)||(rcondR<minrcondR))
        SIGMA=((eps'*eps)^2)/(eps'*R*eps)/Nfree;
        CRIT=SIGMA;
        if(~isfinite(CRIT))
            rcondR
            save DEBUG_BUILD_MINI_EMULATOR2;
            disp('ERROR: CRIT not finite!!!! Workspace has been saved to DEBUG_BUILD_MINI_EMULATOR2.mat\n');
            bob; %cause the code to crash
        end
        return; 
    end    
    
    [CRIT,Grad,Hess,SIGMA]=...
        CritGradHessWRTTheta(Nfree,eps,Z,THETA,R);
    if(~isfinite(CRIT))
        rcondR
        save DEBUG_BUILD_MINI_EMULATOR3;
        disp('ERROR: CRIT not finite!!!! Workspace has been saved to DEBUG_BUILD_MINI_EMULATOR3.mat\n');
        bob; %cause the code to crash
    end
    
    if(SIGMA==0)
        %disp(sprintf('exit %g inewt=%g',2,0));
        return;
    end

    if((SIGMA<0)||...
        ~isempty(find(~isfinite(Hess),1))||...
        ~isempty(find(~isfinite(Grad),1)))
        save DEBUG_BUILD_MINI_EMULATOR4;
        disp('ERROR: negative SIGMA or non finite Hess OR Grad!!!! Workspace has been saved to DEBUG_BUILD_MINI_EMULATOR4.mat\n');
        bob; %cause the code to crash
    end

    for inewt=1:Nnewt
        OldCRIT=CRIT;              

        if(rcond(Hess)<=SmallNearlySingular)
            disp('nearly singular');
            D=-pinv(Hess)*Grad; %to minimize
            %D=pinv(Hess)*Grad; %to maximize
        else
            D=-Hess\Grad; %to minimize
            %D=Hess\Grad; %to maximize
        end
        D=D*sign(D'*Hess*D);
        D(find(THETA<=MinTheta))=0;
        %D(find(THETA>=MaxTheta))=0;
        lambdamax=(D'*D)^0.5;
        if(lambdamax<=0)
            break;
        end
        D=D/lambdamax;
        
        [D,temp]=prevent_rcondR_less_than_minrcondR(...
            Z,minrcondR,R,rcondR,eigvect,eigval,d1rcondRdti,Q1,Q2,D);
            
        if(temp>=0)
            %disp('temp>0');
            if(temp<lambdamax)
                %disp('rcondR limiting step size');
                lambdamax=temp;
            end
        end

        %keep newton step from taking me into region with
        %%THETA>MaxTheta
        %i=find(D>0);
        %if(~isempty(i))
        %    temp=min((MaxTheta(i)-THETA(i))./(D(i)')); %D(i) is positive
        %    %temp=min((MaxTheta   -THETA(i))./(D(i)')); %D(i) is positive
        %    if(temp<lambdamax)
        %        lambdamax=temp; 
        %    end
        %end

        %keep newton step from taking me into region with
        %THETA<MinTheta
        i=find(D<0);
        if(~isempty(i))
            %THETA is positive, D(i) is negative therefore -max of that is
            %positive
            temp=min((MinTheta(i)-THETA(i))./(D(i)'));  %nonuniform MinTheta

            %temp=min((MinTheta   -THETA(i))./(D(i)')); %for same MinTheta
            %for all dimension            
            if(temp<lambdamax)
                lambdamax=temp;
            end
        end
        
        if(((rcondR>=minrcondR)&&(lambdamax==0))||(lambdamax<0))
            break;
        end

        %save DEBUGMEYADA3;
        THETA=THETA+lambdamax*D';
        %save DEBUGMEYADA4;        
        [THETA,R,rcondR,eigvect,eigval,d1rcondRdti,Q1,Q2]=...
            increase_rcondR_to_minrcondR(Z,minrcondR,THETA);
        %save DEBUGMEYADA5;

        %[Rinv,ifsing]=chol(R);
        %if(ifsing)
        %    save DEBUG_BUILD_MINI_EMULATOR5;
        %    disp('ERROR: singular R!!!! Workspace has been saved to DEBUG_BUILD_MINI_EMULATOR5.mat\n');
        %    bob; %cause the code to crash
        %end
        %   
        %Rinv=Rinv\(Rinv'\I);
        
        
        if(any(THETA<MinTheta))
        %if(any(THETA<MinTheta)||any(THETA>MaxTheta))
            break; 
        end    
        
        DeltaCRIT=abs(OldCRIT-CRIT);
        if((CRIT==0)||(DeltaCRIT==0))
            %disp('break 3');
            break;
        end
        %OldCRIT=CRIT;
        
        if(inewt==1)
            FirstDeltaCRIT=DeltaCRIT;
        elseif(DeltaCRIT<FirstDeltaCRIT/1024) %/1048576)
            %disp('converged')
            %disp('break 5');
            break
        end

        if(DeltaCRIT<CRIT/1048576) %/1073741824) %2^30 approx 10^9
            %disp('break 6');
            break;
        end

        if(inewt<Nnewt)
            [CRIT,Grad,Hess,SIGMA]=...
                CritGradHessWRTTheta(Nfree,eps,Z,THETA,R);

            if(~isfinite(CRIT))
                rcondR
                save DEBUG_BUILD_MINI_EMULATOR6;
                disp('ERROR: CRIT not finite!!!! Workspace has been saved to DEBUG_BUILD_MINI_EMULATOR6.mat\n');
                bob; %cause the code to crash
            end

        end
    end %for inewt=1:Nnewt
    
    SIGMA=((eps'*eps)^2)/(eps'*R*eps)/Nfree;
    CRIT=SIGMA;
    
    if(~isfinite(CRIT))
        rcondR
        save DEBUG_BUILD_MINI_EMULATOR7;
        disp('ERROR: CRIT not finite!!!! Workspace has been saved to DEBUG_BUILD_MINI_EMULATOR7.mat\n');
        bob; %cause the code to crash
    end
    
    warning('on','MATLAB:nearlySingularMatrix');
    warning('on','MATLAB:illConditionedMatrix');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Guess a random theta in a neighborhood that is usually "pretty good"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta=rand_theta(CorLenBnds)
  l2=(CorLenBnds.lmin.*2.^(CorLenBnds.Pmax.*rand(1,CorLenBnds.Ndimin))).^2;
  theta=1.0./(2*l2);  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given a theta, compute R, and rcondR (reciprocal of condition number of R)
%if rcondR is below minrcondR, change theta in direction of steepest ascent
%of rcondR, until rcondR >= minrcondR.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [THETA,R,rcondR,eigvect,eigval,d1rcondRdti,Q1,Q2]=...
        increase_rcondR_to_minrcondR(Z,minrcondR,THETA)
    %save DEBUGMEYADA6;
    %disp('hello: increase_rcondR_to_minrcondR')
    %oldtoc=toc;
    Ndimin=size(THETA,2);
    Nsample=sqrt(size(Z,1));
    I=eye(Nsample);
    
    d1rcondRdti=zeros(1,Ndimin);
    Niter=6;
    for iter=1:Niter+1
        %toc1=toc
        R=reshape(exp(Z*THETA'),Nsample,Nsample);
    
        %save DEBUGMEYADA7;
        [eigvect,eigval]=eig(R);
        eigval=diag(eigval);
        [garb,imin]=min(eigval);
        [garb,imax]=max(eigval);
        
        eigvect=eigvect(:,[imin imax]);
        eigval =eigval([imin imax]);
        
        rcondR=eigval(1)/eigval(2);

        if(iter==1)
            OLDrcondR=rcondR;
        end
        
        for idimin=1:Ndimin
            d1Rdti=reshape(Z(:,idimin),Nsample,Nsample).*R;
            d1eigvaldti=[eigvect(:,1)'*d1Rdti*eigvect(:,1) ...
                         eigvect(:,2)'*d1Rdti*eigvect(:,2)];
            d1rcondRdti(idimin)=...
                (d1eigvaldti(1)*eigval(     2)-...
                 eigval(     1)*d1eigvaldti(2) ...
                )/eigval(2)^2;
        end

        0/0;
        warning('off','MATLAB:SingularMatrix');        
        Qtemp=R-eigval(1)*I+eigvect(:,1)*eigvect(:,1)';
        Q1=inv(Qtemp);
        [blah,blahblah]=lastwarn;
        if(numel(findstr(lower(blahblah),'singularmatrix')))
            Q1=pinv(Qtemp);
        end
        0/0;
        Qtemp=R-eigval(2)*I+eigvect(:,2)*eigvect(:,2)';
        Q2=inv(Qtemp);
        [blah,blahblah]=lastwarn;
        if(numel(findstr(lower(blahblah),'singularmatrix')))
            Q2=pinv(Qtemp);
        end
        warning('on','MATLAB:singularMatrix');
        
        %toc5=toc
        if(~(rcondR<minrcondR)||(iter==Niter))
            break;
        end
        
        lambdamax=norm(d1rcondRdti);
        D=d1rcondRdti/lambdamax;
        
        dneigval=eigen_rts_dist(Nsample,Ndimin,I,Z,R,eigval,eigvect,Q1,Q2,D);
        %toc7=toc
        yada=[-1 minrcondR]*dneigval;
        Maxorder=length(dneigval)-1;
        for maxorder=Maxorder:-1:1
            pf=factorial(maxorder:-1:0);
            %save DEBUGMEWTF;
            rts=roots(yada./pf);
            
            realrts=real(rts(:));
            i=find(isreal(rts)&(realrts>=0));
            
            if(~isempty(i))
                lambdamax=1.05*min(realrts(i));
                break;
            end
            yada=yada(2:length(yada));
        end
        %toc10=toc  
        i=find(D<0);
        if(~isempty(i))
            lambdamax=min(0.95*min((-THETA(i))./(D(i))),lambdamax);
        end

        lambdamax=max(lambdamax,0);
        if(lambdamax==0)
            break;
        end
        
        %save DEBUGTHISQUICK;
        THETA=THETA+D*lambdamax;
        
        %THETA=THETA+(minrcondR-rcondR)/(d1rcondRdti*D')*D;
    end
    %disp(sprintf('goodbye: increase_rcondR_to_minrcondR: %g [sec]',toc-oldtoc));
    %bob
    %disp(sprintf('iter=%g rcondR/minrcondR=%g/%g=%g OLDwas=%g/%g=%g',iter,...
    %    rcondR,minrcondR,rcondR/minrcondR,OLDrcondR,minrcondR,OLDrcondR/minrcondR));
    %bob;
return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the [4th 3rd 2nd 1st 0th] derivatives (with respect to distance
%traveled in the direction of vector D) of the largest and smallest
%eigenvalues of R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The implementation assumes R is real symmetric postitive definite the math
%comes from...
%
%Jankovic, Exact nth Derivatives of Eigenvalues and Eigenvectors.
%Journal of Guidance, Control, and Dynamics, Vol 17. No 1,
%January-February 1994 pp 136-144
%
%Note that I don't normalize, i.e. normalize by the identity matrix => K=I
%my matrix R is real symmetric and positive definite => complex conjugate 
%of an eigenvector is the eigenvector itself.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dneigval=eigen_rts_dist(Nsample,Ndimin,I,Z,R,eigval,eigvect,Q1,Q2,D)

    d1R=zeros(Nsample);
    d2R=d1R;
    d3R=d1R;
%     d4R=d1R;
    %tocA=toc
    for idimin=1:Ndimin
        temp1=reshape(Z(:,idimin),[Nsample Nsample]).*R*D(idimin);
        d1R=d1R+temp1;
        for jdimin=1:Ndimin
            temp2=reshape(Z(:,jdimin),[Nsample Nsample]).*temp1*D(jdimin);
            d2R=d2R+temp2;
%             for kdimin=1:Ndimin
%                 temp3=reshape(Z(:,kdimin),[Nsample Nsample]).*temp2*D(kdimin);
%                 d3R=d3R+temp3;
% %                 for ldimin=1:Ndimin
% %                     temp4=reshape(Z(:,ldimin),[Nsample Nsample]).*temp3*D(ldimin);
% %                     d4R=d4R+temp4;
% %                 end
%             end
        end
    end
    %tocB=toc
    z1       =d1R*eigvect;
    d1eigval =[eigvect(:,1)'*z1(:,1) eigvect(:,2)'*z1(:,2)];
    eta1     =[0                     0                    ];
    d1eigvect=[...
        -0.5*eta1(1)*eigvect(:,1)-Q1*(eigvect(:,1)*d1eigval(1)+z1(:,1))...
        -0.5*eta1(2)*eigvect(:,2)-Q2*(eigvect(:,2)*d1eigval(2)+z1(:,2))];
    
    z2       =d2R*eigvect+[...
        2*(d1R-d1eigval(1)*I)*d1eigvect(:,1)...
        2*(d1R-d1eigval(2)*I)*d1eigvect(:,2)];
    d2eigval =[eigvect(:,1)'*z2(:,1) eigvect(:,2)'*z2(:,2)];
    eta2     =2*sum(d1eigvect.^2);
    d2eigvect=[...
        -0.5*eta2(1)*eigvect(:,1)-Q1*(eigvect(:,1)*d2eigval(1)+z2(:,1))...
        -0.5*eta2(2)*eigvect(:,2)-Q2*(eigvect(:,2)*d2eigval(2)+z2(:,2))];
    
%     z3       =d3R*eigvect+[...
%         3*(d2R-d2eigval(1)*I)*d1eigvect(:,1)+3*(d1R-d1eigval(1)*I)*d2eigvect(:,1)...
%         3*(d2R-d2eigval(2)*I)*d1eigvect(:,2)+3*(d1R-d1eigval(2)*I)*d2eigvect(:,2)];
%     d3eigval =[eigvect(:,1)'*z3(:,1) eigvect(:,2)'*z3(:,2)];
%     eta3     =6*sum(d2eigvect.*d1eigvect);
%     d3eigvect=[...
%         -0.5*eta3(1)*eigvect(:,1)-Q1*(eigvect(:,1)*d3eigval(1)+z3(:,1))...
%         -0.5*eta3(2)*eigvect(:,2)-Q2*(eigvect(:,2)*d3eigval(2)+z3(:,2))];
    
%     z4       =d4R*eigvect+[...
%         4*(d3R-d3eigval(1)*I)*d1eigvect(:,1)+6*(d2R-d2eigval(1)*I)*d2eigvect(:,1)+4*(d1R-d1eigval(1)*I)*d3eigvect(:,1)...
%         4*(d3R-d3eigval(2)*I)*d1eigvect(:,2)+6*(d2R-d2eigval(2)*I)*d2eigvect(:,2)+4*(d1R-d1eigval(2)*I)*d3eigvect(:,2)];
%     d4eigval =[eigvect(:,1)'*z4(:,1) eigvect(:,2)'*z4(:,2)];
    %tocG=toc
%    dneigval=[d4eigval(:) d3eigval(:) d2eigval(:) d1eigval(:) eigval(:)];
%    dneigval=[d3eigval(:) d2eigval(:) d1eigval(:) eigval(:)];
    dneigval=[d2eigval(:) d1eigval(:) eigval(:)];
    %save DEBUGMEWTF2;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make sure the step to minimize the criteria doesn't take rcondR below
%minrcondR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,lambdamax]=prevent_rcondR_less_than_minrcondR(...
    Z,minrcondR,R,rcondR,eigvect,eigval,d1rcondRdti,Q1,Q2,D)
 
    %if(rcondR<minrcondR)
    %    lambdamax=0;
    %    return;
    %end
    %normdir=d1rcondRdti'/norm(d1rcondRdti); 

    %dotproduct=D'*normdir;
    %if((rcondR<=1.05*minrcondR)&&(dotproduct<0))
    %    D=D-dotproduct*normdir;
    %end

    Nsample=size(R,1);
    Ndimin=size(Z,2);
   
    I=eye(Nsample);

    dneigval=eigen_rts_dist(Nsample,Ndimin,I,Z,R,eigval,eigvect,Q1,Q2,D);
    yada=[-1 minrcondR]*dneigval;
    maxorder=length(dneigval)-1;
    pf=factorial(maxorder:-1:0);

    %save DEBUGMEWTF;
        
    rts=roots(yada./pf);


    
    realrts=real(rts);
    i=find(realrts>=0);
    if(~isempty(i))
        lambdamax=0.95*min(realrts(i));
    elseif(rcondR<minrcondR)
        lambdamax=0;
    else
        lambdamax=inf;
        %
        %save DEBUG_BUILD_MINI_EMULATOR8;
        %disp('ERROR: rcond>=minrcondR and no non-negative lambdamax!!!! Workspace has been saved to DEBUG_BUILD_MINI_EMULATOR8.mat\n');
        %bob; %cause the code to crash
    end
   
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the "criteria," its gradient and Hessian, and Sigma, just in case
%the criteria is different than Sigma (unadjusted variance), this will be
%used in the next Newton Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I've stripped out Busby's calcuation of Sigma (unadjusted variance)
%Right now the criteria to minimize is sigma, might want try minimizing the
%sum over 1 hop simplices of adjusted variance at the simplex centroid
%weighted by the simplex "content" (hypervolume)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CRIT,Grad,Hess,Sigma]=...
    CritGradHessWRTTheta(Nfree,eps,Z,theta,R)

    Nsample=size(eps,1);
    Ndimin=size(theta,2);
    %Nbasis=size(G,2), if I need Nbasis I'll have to pass it in;

    d1Sigmafast=zeros(Ndimin,1);
    d2Sigmafast=zeros(Ndimin);

    %Rinv_eps=Rinv*eps;
    epst_eps=eps'*eps;
    if(epst_eps^4==0)
        Sigmafast=0;
        CRIT=Sigmafast;
        Sigma=Sigmafast;
        Grad=d1Sigmafast;
        Hess=d2Sigmafast;
        return;
    end

    Sigmafast=((epst_eps)^2)/(eps'*R*eps);  %actuall Sigmafast*Nsample

    %d1R_Rinv_eps=zeros(Nsample,Ndimin);
    for idimin=1:Ndimin
        dRdti=reshape(Z(:,idimin),[Nsample Nsample]).*R;
        
        d1Sigmafast(idimin)=...
            -(Sigmafast/epst_eps)^2*(eps'*dRdti*eps);
   
        %d1R_Rinv_eps(:,idimin)=dRdti*Rinv_eps;
    end

    temp=reshape(Z(:,1),[Nsample Nsample]);
    dRdti=temp.*R;
    dRdtj=dRdti;
    istart=1;
    for idimin=1:Ndimin
        jdimin=idimin; %could instead say jdimin=1 before this loop 
        d2Rdtidtj=temp.*dRdti;  

        d2Sigmafast(idimin,jdimin)=...
            +2*Sigmafast^3/epst_eps^4*(eps'*dRdtj*eps)*(eps'*dRdti*eps)...
            -(Sigmafast/epst_eps)^2*(eps'*d2Rdtidtj*eps);

        for jdimin=Ndimin:-1:idimin+1
            temp=reshape(Z(:,jdimin),[Nsample Nsample]);
            dRdtj=temp.*R;
            d2Rdtidtj=temp.*dRdti;
            
            d2Sigmafast(idimin,jdimin)=...
                +2*Sigmafast^3/epst_eps^4*(eps'*dRdtj*eps)*(eps'*dRdti*eps)...
                -(Sigmafast/epst_eps)^2*(eps'*d2Rdtidtj*eps);

            d2Sigmafast(jdimin,idimin)=d2Sigmafast(idimin,jdimin);
        end
        dRdti=dRdtj;
        
        istart=istart+Ndimin-(idimin-1);
    end

      Sigmafast=  Sigmafast/Nfree;
    d1Sigmafast=d1Sigmafast/Nfree;
    d2Sigmafast=d2Sigmafast/Nfree;

    Sigma=Sigmafast;
    CRIT=Sigmafast;
    Grad=d1Sigmafast;
    Hess=d2Sigmafast;
return;
