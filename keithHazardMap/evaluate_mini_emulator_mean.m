function evaluate_mini_emulator_mean(samplenumber,macrosimplex_resample_filename)
    iftoc=1;
    if(iftoc)
        tic;
    end

    if(ischar(samplenumber))
        samplenumber=str2double(samplenumber);
    end

    Nmoments=1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read input file #1: this contains the mini-emulator
    miniemulatorfilename=sprintf('mini_emulator.%06d',samplenumber);
    fid=fopen(miniemulatorfilename,'r');
    
    Nskip=sscanf(fgets(fid),'additional file format lines=%g',1);
    for i=1:Nskip
        fgets(fid);
    end

%         %this data can be printed here
%         fprintf(fid,'%%iymacro: index/id of the simulation at the center of this mini-emulator, this should be the same as the filename suffix: [1] integer\n');
%         fprintf(fid,'%%Ndiminspatial=2: #number of spatial input dimension (east,north): [1] integer\n');
%         fprintf(fid,'%%Ndimin=6: # of input dimensions (east,north,log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
%         fprintf(fid,'%%Ny0T: #of data points (y0) in this mini-emulator''s central simulation: [1] integer\n');
%         fprintf(fid,'%%Ny0N: #of micro-emulators in this mini-emulator: [1] integer\n');
%         fprintf(fid,'%%Ny: #of data points in this mini-emulator: [1] integer\n');        
%         fprintf(fid,'%%{{y=(east,north,log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg])},{sy=h}} [Ndimin+1] doubles, there are [Ny] of these lines, the first [Ny0N] of them are the centers of the micro-emulators\n');
%         fprintf(fid,'%%NsimpT: # of simplices in y0 tesselation: [1] integer\n');
%         fprintf(fid,'%%NsimpN: # of simplices in y0 tesselation that have at least 1 micro-emulator as a node: [1] integer\n');
%         fprintf(fid,'%%tess0: [NsimpT Ndimin+1] array of integers/data-point-indices\n');
%         
%         %this data must be printed in the following loop
%         fprintf(fid,'%%then there will be [Ny0N] micro-emulator "structures," each of which contains the following...\n');
%         fprintf(fid,'%%(iy0): index of micro emulator, should also be the first entry in iy0neigh and iyneigh: [1] integer\n');
%         fprintf(fid,'%%lmax: maximum plausible correlation lengths: [Ndimin] doubles\n');
%         fprintf(fid,'%%lmin: minimum plausible correlation lengths: [Ndimin] doubles\n');
%         fprintf(fid,'%%theta: roughness parameters: [Ndimin] doubles\n');
%         fprintf(fid,'%%sigma2: unadjusted variance: [1] double\n');
%         fprintf(fid,'%%beta: least squares coefficients: [Ndimin+1] doubles\n');
%         fprintf(fid,'%%Nneigh0: #of y0 data points in the neighborhood of this micro-emulator: [1] integer\n');
%         fprintf(fid,'%%{{iy0neigh: index of y0 neighbor: [1] integer},{i0hop: #of y0 hops from center of micro-emulator: [1] integer}}\n');
%         fprintf(fid,'%%Nneigh: #of y data points in the neighborhood of this micro-emulator: [1] integer\n');
%         fprintf(fid,'%%{{iyneigh: index of y neighbor: [1] integer},{ihop: #of (y0,ymacro) hops from center of micro-emulator: [3] integers},{Rinveps: R^-1*eps (eps = sy minus unadjusted micro-emulator mean): [1] double}}\n');
    

    checksamplenumber=sscanf(fgets(fid),'%g',1);
    if(checksamplenumber~=samplenumber)
        save DEBUG_EVALUATE_MINI_EMULATOR_MEAN1;
        disp(sprintf('ERROR: file "%s" claims to be mini-emulator #%g instead!!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN1.mat\n',miniemulatorfilename,checksamplenumber));
        bob; %cause the code to crash
    end

    Ndiminspatial=sscanf(fgets(fid),'%g',1);
    Ndimin=sscanf(fgets(fid),'%g',1);
    Ny0T=sscanf(fgets(fid),'%g',1);
    Ny0N=sscanf(fgets(fid),'%g',1);
    Ny=sscanf(fgets(fid),'%g',1);
    
    yada=sprintf('%%g%s\n',repmat(' %g',1,Ndimin)); %want it to crash if the data isn't in this EXACT format
    ysy=fscanf(fid,yada,[Ndimin+1 Ny])';

    %ysy(Ny,:)
    
    %fgets(fid); %read in the \n that fscanf left on the end of the
    %previous line: DANGER MATLAB IS NOT PERFORMING AS EXPECTED, it isn't
    %leaving the \n on the previous line
    
    Ntess0T=sscanf(fgets(fid),'%g',1);
    Ntess0N=sscanf(fgets(fid),'%g',1);
    yada=sprintf('%%g%s\n',repmat(' %g',1,Ndiminspatial)); %want it to crash if the data isn't in this EXACT format
    tess0=fscanf(fid,yada,[Ndiminspatial+1 Ntess0T])';

    %sizetess0=size(tess0)
    
    %fgets(fid); %read in the \n that fscanf left on the end of the previous line    
    
    thetaformat=sprintf('%%g%s',repmat(' %g',1,Ndimin-1));
    betaformat=sprintf('%%g%s',repmat(' %g',1,Ndimin));
    
    %only allocate space for the parts of the micro-emulators needed to
    %evaluate the adjusted mean (should do a similar minimization of 
    %retained data for evaluating the adjusted variance)
    pwem=repmat(struct('iyneigh',[],'thetaTRAN',[],'beta',[],'Rinveps',[]),Ny0N,1);
            
    for iy0=1:Ny0N
        %blah=fgets(fid)
        %checkiy0=sscanf(blah,'(%g)',1)
        checkiy0=sscanf(fgets(fid),'(%g)',1);
        if(isempty(checkiy0)||(iy0~=checkiy0))
            save DEBUG_EVALUATE_MINI_EMULATOR_MEAN2;
            disp(sprintf('ERROR: iy0=%g ~= checkiy0=%g!!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN2.mat\n',iy0,checkiy0));
            bob; %cause the code to crash
        end

        lmax=sscanf(fgets(fid),thetaformat,Ndimin)'; %don't need to retain this for adjusted mean
        lmin=sscanf(fgets(fid),thetaformat,Ndimin)'; %don't need to retain this for adjusted mean
        
        pwem(iy0).thetaTRAN=sscanf(fgets(fid),thetaformat,Ndimin);        
        sigma=sscanf(fgets(fid),'%g',1); %don't need to retain this for adjusted mean
        pwem(iy0).beta =sscanf(fgets(fid),betaformat,Ndimin+1);
        
        Nneigh0=sscanf(fgets(fid),'%g',1); %don't need to retain this for adjusted mean
        yada=fscanf(fid,'%g %g\n',[2 Nneigh0])'; %[iy0neigh i0hop], don't need to retain these for adjusted mean

        %fgets(fid) %read in the \n that fscanf left on the end of the previous line    
            
        Nneigh=sscanf(fgets(fid),'%g',1);  %don't need to retain this for adjusted mean
        yada=fscanf(fid,'%g %g %g %g %g\n',[5 Nneigh])'; %[iyneigh ihop Rinveps], don't need to retain ihop for adjusted mean
        pwem(iy0).iyneigh=yada(:,1);
        pwem(iy0).Rinveps=yada(:,5);
        %fgets(fid); %read in the \n that fscanf left on the end of the previous line  
        %bob
    end
    fclose(fid);
    clear lmax lmin sigma Nneigh0 Nneigh yada;
    
    if(iftoc)
        disp(sprintf('Done reading in mini-emulator at time t=%g [sec]',toc));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read input file #2: this contains points at which to evaluate the
    %mini-emulator

    fid=fopen(macrosimplex_resample_filename,'r');

    Nskip=sscanf(fgets(fid),'additional file format lines=%g',1);
    for i=1:Nskip
        fgets(fid);
    end

%         fprintf(fid,'%%Unique Key: the number before the decimal is the simplex id. The number AFTER THE DECIMAL is randomly generated after seeding by clock and also occurs in mini-emulator evaluation output filenames allowing you to match them: [1] double\n');
%         fprintf(fid,'%%Ndiminmacro=4: #of uncertain dimensions, (log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]) for each resample input: [1] integer\n');
%         fprintf(fid,'%%{{(simplex id): a sanity check: [1] integer},{id''s of simulation/mini-emulator in the simplex: [Ndiminmacro+1] integers}}\n');
%         fprintf(fid,'%%Ndiminspatial=2: #of spatial dimensions (east,north): [1] integer\n');
%         fprintf(fid,'%%Nmacroseachx: #of macro input points associated with each spatial point, done to reduce file size but remain general: [1] integer\n');
%         fprintf(fid,'%%Nxspatial: #of LISTED spatial coordinates: [1] integer\n');
%         fprintf(fid,'%%file format option indicator: [1] integer\n');
%         fprintf(fid,'%%one of the following two options\n');
%         fprintf(fid,'%%option 1: for when all spatial points need to be evaluated at the same set(s) of macro coordinates\n');
%         fprintf(fid,'%%     one line containing: [Nmacroseachx] sets of {(log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): macro input coordinates: [Ndiminmacro] doubles}\n');
%         fprintf(fid,'%%     [Nxspatial] lines containing {(east,north): [2] doubles}\n');
%         fprintf(fid,'%%option 2: for when each spatial point can be paired with different macro input(s), in this case there can also be duplicates of listed (east,north) coordinates that have different macro coordinates\n');
%         fprintf(fid,'%%     [Nxspatial] lines containing {{(east,north): [2] doubles},{[Nmacroseachx] sets of {(log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): macro input coordinates: [Ndiminmacro] doubles}}}\n');
    
    yada=sscanf(fgets(fid),'%d.%d',2);
    iMacroSimplex=yada(1);
    RandomKey=yada(2);
    UniqueMacroSimplexKey=sprintf('%d.%08d',iMacroSimplex,RandomKey);
    UniqueMiniEmulatorKey=sprintf('%s.1.%06d',UniqueMacroSimplexKey,samplenumber); %the 1 in the .1. indicates it's an evaluation of the adjusted first moment (mean)
    
    Ndiminmacro=sscanf(fgets(fid),'%g',1);
    if(Ndimin-(Ndiminspatial+Ndiminmacro))
        save DEBUG_EVALUATE_MINI_EMULATOR_MEAN3;
        disp('ERROR: Ndimin~=Ndiminspatial+Ndiminmacro!!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN3.mat\n');
        bob; %cause the code to crash
    end
    
    yada=sscanf(fgets(fid),sprintf('(%%g)%s',repmat(' %g',1,Ndiminmacro+1)),Ndiminmacro+2);
    if(yada(1)~=iMacroSimplex)
        save DEBUG_EVALUATE_MINI_EMULATOR_MEAN4;
        disp(sprintf('ERROR: file "%s" claims it''s for simplex %g!!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN4.mat\n',macrosimplex_resample_filename,yada(1)));
        bob; %cause the code to crash
    end
    
    MacroSimplexNodes=yada(2:Ndiminmacro+2);
    iMacroSimplexNode=find(MacroSimplexNodes==samplenumber);
    if(isempty(iMacroSimplexNode))
        save DEBUG_EVALUATE_MINI_EMULATOR_MEAN5;
        disp(sprintf('ERROR: file "%s" claims its only macro simplex nodes are %s but this mini-emulator (#%g) is not one of them!!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN5.mat\n',...
            macrosimplex_resample_filename,sprintf(sprintf('{%%g%s}',repmat(',%g',Ndiminmacro)),MacroSimplexNodes),samplenumber));
        bob; %cause the code to crash
    end
       
    checkNdiminspatial=sscanf(fgets(fid),'%g',1);
    if(Ndiminspatial~=checkNdiminspatial)
        save DEBUG_EVALUATE_MINI_EMULATOR_MEAN6;
        disp(sprintf('ERROR: file Ndiminspatial=%g ~= checkNdiminspatial=%g!!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN6.mat\n',...
            Ndiminspatial,checkNdiminspatial));
        bob; %cause the code to crash
    end
    
    Nmacroseachx=sscanf(fgets(fid),'%g',1);
    if(~(Nmacroseachx>0))
        save DEBUG_EVALUATE_MINI_EMULATOR_MEAN7;
        disp(sprintf('ERROR: Nmacroseachx=%g is not greater than zero!!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN7.mat\n',Nmacroseachx));
        bob; %cause the code to crash
    end
            
    Nxspatial=sscanf(fgets(fid),'%g',1);
    fileformatoption=sscanf(fgets(fid),'%g',1);

    Iy0T=1:Ny0T;
    
    switch(fileformatoption)
        case 1,
            Xmacro=sscanf(fgets(fid),'%g',[Ndiminmacro Nmacroseachx])';
            yada=sprintf('%%g%s\n',repmat(' %g',1,Ndiminspatial-1));
            xspatial=fscanf(fid,yada,[Ndiminspatial Nxspatial])';
            T=repmat(tsearch(ysy(Iy0T,1),ysy(Iy0T,2),tess0,xspatial(:,1),xspatial(:,2)),Nmacroseachx,1);        
            Nx=Nxspatial*Nmacroseachx;
            %save DEBUGMEWTFevalmini
            xsx=[repmat(xspatial,Nmacroseachx,1) reshape(repmat(Xmacro',Nxspatial,1),Ndiminmacro,Nxspatial*Nmacroseachx)' zeros(Nx,1)];
            clear xspatial Nxspatial Nmacroseachx;
        case 2,            
            yada=sprintf('%%g%s\n',repmat(' %g',1,Ndiminspatial-1+Nmacroseachx*Ndiminmacro));
            yada=fscanf(fid,yada,[Ndiminspatial+Nmacroseachx*Ndiminmacro Nxspatial])';

            if(Nmacroseachx==1)
                Nx=Nxspatial;
                xsx=[yada zeros(Nx,1)];
                clear yada Nxspatial Nmacroseachx;
                T=tsearch(ysy(Iy0T,1),ysy(Iy0T,2),tess0,xsx(:,1),xsx(:,2));
            else
                xspatial=yada(:,1:Ndiminspatial);
                T=repmat(tsearch(ysy(Iy0T,1),ysy(Iy0T,2),tess0,xspatial(:,1),xspatial(:,2)),Nmacroseachx,1);
                yada=yada(:,Ndiminspatial+(1:Nmacroseachx*Ndiminmacro));
                ii=repmat((0:Nmacroseachx-1)*Ndiminmacro,1,Ndiminmacro)+...
                    reshape(repmat(1:Ndiminmacro,Nmacroseachx,1),1,Ndiminmacro*Nmacroseachx);
                Nx=Nxspatial*Nmacroseachx;
                xsx=[repmat(xspatial,Nmacroseachx,1) ...
                     reshape(yada(:,ii),Nx,Ndiminmacro) zeros(Nx,1)];                
                clear yada xspatial Nxspatial Nmacroseachx ii;
            end
        otherwise,
            save DEBUG_EVALUATE_MINI_EMULATOR_MEAN8;
            disp(sprintf('ERROR: unknown 2nd input file format option!!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN8.\n'));
            bob; %cause the code to crash
    end   

    
    if(iftoc)
        disp(sprintf('Done reading in resamples at time t=%g [sec]',toc));
    end
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now evaluate the min-emulator's adjusted mean 

    Idimout=Ndimin+1;
    Nnodes=Ndiminspatial+1;
    
    [T,isort]=sort(T);    
    if(iftoc)
        disp(sprintf('done sorting T at time t=%g [sec]',toc));
    end
            
    xsx=xsx(isort,:);
    
    ilastinflow=find(T<=Ntess0N,1,'last');
        
    Iend=find(diff(T(1:ilastinflow)));
    Istart=[1; Iend+1];
    Iend=[Iend; ilastinflow];
    Nmicrosimptoeval=numel(Istart);
    ONES3_1=ones(3,1);
    for imicrosimptoeval=1:Nmicrosimptoeval
        istart=Istart(imicrosimptoeval);
        iend  =Iend(  imicrosimptoeval);
        ixS=istart:iend; %indices of x in this micro Simplex
        isimpmicro=T(istart);
        if(iftoc)
            disp(sprintf('isimpmicro=%g at time t=%g [sec]',isimpmicro,toc));
        end
        Iy0nodes=tess0(isimpmicro,:);
        
        xsxS=xsx(ixS,:);
        P=[xsxS(:,1:2) ones(numel(ixS),1)]/[ysy(Iy0nodes,1:2) ONES3_1]; %new code so I can use tsearch instead of tsearchn
        

        for inode=1:Nnodes
            ipwem=Iy0nodes(inode);
             xsxS=add_micro_emulator_adjusted_mean(pwem(ipwem),...
                ysy(pwem(ipwem).iyneigh,:),xsxS,P(:,inode));
%                ysy(pwem(ipwem).iyneigh,:),xsxS,P(ixS,inode)); %for use
%                with tsearchn
        end
        xsx(ixS,Idimout)=xsxS(:,Idimout);
    end
    
    if(iftoc)
        disp(sprintf('Done evaluating resamples at time t=%g [sec]',toc));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now write the min-emulator's output (adjusted mean) to the datafile

    xsx=sortrows(xsx,[Ndimin:-1:1 Idimout]);    
    idiffmacro=[find(any(diff(xsx(:,Ndiminspatial+1:Ndimin),1,1),2)); Nx];
    Nmacroseachx=numel(idiffmacro);
    Nxmap=unique([idiffmacro(1); diff(idiffmacro)]);
    
    if((numel(Nxmap)==1)&&(Nmacroseachx>1))    
        fileformatoption=1;
    else
        fileformatoption=2;
    end
    
    filenameout=sprintf('mini_emulator_eval_%s',UniqueMiniEmulatorKey);
    fid=fopen(filenameout,'w');
    
    Nskip=16;
    fprintf(fid,'additional file format lines=%g\n',Nskip);
    
    fprintf(fid,'%%a UniqueKey that contains simplex id, the randomly generated after seeded by time key passed in (allows for pairing), the degree of the statistical moment evaluated (1=mean, 2=variance, 3=skewness, 4=kurtosis,...,-2=first 2 moments, -3=first 3 moments), and the 6 digit (left padded by zeros) id of this mini-emulator\n');
    fprintf(fid,'%%Nmoments: Nmoments=1 if the moment key is positive, otherwise it is absolute value of the moment key: [1] integer\n');
    fprintf(fid,'%%Ndiminmacro=2: #of uncertain dimensions (log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]), for each resample input: [1] integer\n');
    fprintf(fid,'%%{{(simplex id): a sanity check: [1] integer},{id''s of simulation/mini-emulator in the simplex, another sanity check this mini-emulator should one of the ones listed here: [Ndiminmacro+1] integers}}\n');
    fprintf(fid,'%%Ndimin=6: #of mini-emulator input dimensions (east,north,log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
    fprintf(fid,'%%file format option indicator: [1] integer\n');
    fprintf(fid,'%%one of the following two options\n');
    fprintf(fid,'%%option 1: for when all map points need to be evaluated at the same set(s) of macro coordinates\n');
    fprintf(fid,'%%     Nmacroseachx: #of macro input points associated with each map point: [1] integer\n');
    fprintf(fid,'%%     one line containing: [Nmacroseachx] sets of {(log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): macro input coordinates: [Ndiminmacro] doubles}\n');
    fprintf(fid,'%%     Nxmap: #of unique map-points: [1]\n');
    fprintf(fid,'%%     [Nxmap] lines containing {(east,north): [2] doubles}\n');
    fprintf(fid,'%%     [(Nxmap)*Nmoments] lines containing: {adjusted moment for Nmacroseachx macro inputs: [Nmacroseachx] doubles}, if Nmoments >1, all mean values come before all variance values, etc\n'); 
    fprintf(fid,'%%option 2: for when each map point can be paired with different macro input(s), in this case there can also be duplicates of listed (east,north) coordinates that have different macro coordinates\n');
    fprintf(fid,'%%     Nx: #of (Ndimin-tuple) input points the mini-emulator was evalutated at: [1] integer\n');
    fprintf(fid,'%%     [Nx] lines, each of which contains {{x=(east,north,log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [Ndimin] doubles},{the adjusted statistical moment(s): [Nmoment] double(s)}}\n');

    fprintf(fid,'%s\n',UniqueMiniEmulatorKey);
    fprintf(fid,'%d\n',Nmoments,Ndiminmacro);
    fprintf(fid,sprintf('(%d)%s\n',iMacroSimplex,repmat(' %d',1,Ndiminmacro+1)),MacroSimplexNodes);
    fprintf(fid,'%d\n',Ndimin,fileformatoption);
    switch(fileformatoption)
        case 1,
            xmap=xsx(1:Nxmap,1:2);
            macroseachx=xsx(idiffmacro,3:Ndimin);            
            xsx=permute(reshape(xsx(:,Idimout),[Nxmap Nmacroseachx Nmoments]),[2 1 3]);
            fprintf(fid,'%d\n',Nmacroseachx);
            yada=sprintf('%%.10g%s\n',repmat(' %.10g',1,Nmacroseachx*Ndiminmacro-1));
            fprintf(fid,yada,macroseachx');
            fprintf(fid,'%g\n',Nxmap);
            fprintf(fid,'%.2f %.2f\n',xmap');
            yada=sprintf('%%.10g%s\n',repmat(' %.10g',1,Nmacroseachx-1));
            fprintf(fid,yada,xsx);            
        case 2,
            fprintf(fid,'%g\n',Nx);
            yada=sprintf('%%.2f%s%s\n',repmat(' %.2f',1,Ndiminspatial-1),repmat(' %.10g',1,Ndiminmacro+Nmoments));
            fprintf(fid,yada,xsx');            
        otherwise,
            save DEBUG_EVALUATE_MINI_EMULATOR_MEAN9;
            disp(sprintf('ERROR: I don''t know how to handle file format %d!!! Workspace saved to DEBUG_EVALUATE_MINI_EMULATOR_MEAN9.mat\n',fileformatoption));
            bob; %cause the code to crash
    end
    
    fclose(fid);
    
    if(iftoc)
        disp(sprintf('Done writing evaluated resamples at time t=%g [sec]',toc));
    end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%evaluates the Bayes Linear (micro) Emulator AT a collection of points x in
%input space, the UNADJUSTED mean is
% max([1 x(*,1)-y(1,1) x(*,2)-y(1,2) ... x(*,Ndimin)-y(1,Ndimin)]*beta,0)
%
%pwem is the piecewise micro emulator
%
%ysy are the input dimensions and output dimension (singular) of the known 
%data points
%
%xsx are the input dimensions and output dimension (singular) of the points
%to be evaluated.
%
%p is a vector of weights (the barycentric coordinate of x for THIS node of
%the mini-emulator simplex) to multiply the micro-emulator's adjusted mean
%by before adding it to the previous contents of the output, sx, at each of
%the evaluation points, x.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xsx=add_micro_emulator_adjusted_mean(pwem,ysy,xsx,p)
    %if(~any(pwem.beta))
    %    return;
    %end
    
    Idimout=size(ysy,2);
    Ndimin=Idimout-1;
    Idimin=1:Ndimin;
    Ny=size(ysy,1);
    Nx=size(xsx,1);
    NxNy=Nx*Ny;
    ONES=ones(Nx,1);
        
    if(nargin<4)
        %if the user doesn't give you p assume that p is one (he's just
        %evaluating one micro-emulator rather than the all of the micro
        %emulator nodes of a mini-emulator spatial simplex)
        p=ONES;
    end
    
    negZ=zeros(NxNy,Ndimin);
    for idimin=1:Ndimin
       negZ(:,idimin)=reshape((repmat(xsx(:,idimin),1,Ny)-repmat(ysy(:,idimin)',Nx,1)).^2,NxNy,1);
    end

    %if(size(pwem.beta,1)~=Idimout)
    %    save DEBUGME_EVAL_MINI_MEAN_WTF1;
    %bob;
    %end
    
    xsx(:,Idimout)=xsx(:,Idimout)+p.*(...
        max([ONES xsx(:,Idimin)-repmat(ysy(1,Idimin),Nx,1)]*pwem.beta,0)+...
        reshape(exp(negZ*(-pwem.thetaTRAN)),Nx,Ny)*pwem.Rinveps);
return;