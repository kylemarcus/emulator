function extract_mini_emulator_build_meta_data(samplenumber)
    ifplot=0; %1;
    if(ischar(samplenumber))
      samplenumber=str2num(samplenumber);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read the macro emulator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid=fopen('macro_emulator.pwem','r');
    Nskip=sscanf(fgets(fid),'additional file format lines=%g',1);
    for i=1:Nskip
        fgets(fid);
    end
    NdiminmacroY=sscanf(fgets(fid),'%g',1);
    Nymacro=sscanf(fgets(fid),'%g',1);

    pwem=repmat(struct('Nneigh',0,'iyneigh',[],'ihop',[],'Nsimpneigh',0,...
        'isimpneigh',[],'lmin',[]),Nymacro,1);
    
    y=zeros(Nymacro,NdiminmacroY);
    for iy=1:Nymacro
        thisline=fgets(fid);
        yada=sscanf(thisline,sprintf('(%d) %s',iy,...
            repmat(' %g',1,NdiminmacroY+1)),NdiminmacroY+1);
        y(iy,:)=yada(1:NdiminmacroY);
        pwem(iy).isimpneigh=sscanf(thisline,sprintf('(%d)%s%s',iy,...
            sprintf(' %.10g',yada),repmat(' %g',1,yada(NdiminmacroY+1))),...
            yada(NdiminmacroY+1));
        %size(pwem(iy).isimpneigh)
        pwem(iy).Nsimpneigh=yada(NdiminmacroY+1);
        %save DEBUGME
        %bob
    end

    Nsimpmacro=sscanf(fgets(fid),'%g',1);
    yada=sprintf('(%%g)%s\n',repmat(' %g',1,NdiminmacroY+1));
    tess=fscanf(fid,yada,[NdiminmacroY+2 Nsimpmacro]);
    
    if(any(tess(1,:)~=1:Nsimpmacro))
        %need to display an error message
        return;
    end
    tess=tess(2:NdiminmacroY+2,:)';
    fclose(fid);
    
    if(0&&ifplot)
        trimesh(tess,y(:,1),y(:,2));
        title('read in macro emulator','fontsize',14);
    end

    %use the average distance between points as "equivalent distances" the 
    %scale factor is one/equivalent-distance
    
    scalefactor=1.0./((max(y)-min(y))/Nymacro^(1/NdiminmacroY));
    SFM=repmat(scalefactor,NdiminmacroY+1,1);
    Sn=hypersphere_Sn(NdiminmacroY);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %this information is useful for building mini-emulators, not so much 
    %for identifying mini-emulators to evaluate
    %for iy=1:Nymacro
    for iy=samplenumber
       %only for diagnostic, to test if beyond 2 hops works
       %pwem(iy).Nsimpneigh=pwem(iy).Nsimpneigh+1;
       %pwem(iy).isimpneigh(pwem(iy).Nsimpneigh)=1;
       
       yada=tess(pwem(iy).isimpneigh,:); %real
       included=zeros(1,pwem(iy).Nsimpneigh);
       
       iyneigh=iy;
       iynexthop=[];
       ihop=0;
       Nhop=0;
       
       Nleft=sum(~included);
       while(Nleft)
           Nhop=Nhop+1;
           for isimp=1:pwem(iy).Nsimpneigh
               if(~included(isimp))
                   %save DEBUGME;
                   isect=setdiff(yada(isimp,:),iyneigh);
                   lenisect=length(isect);
                   if((0<lenisect)&&(lenisect<=NdiminmacroY))
                       iynexthop=[iynexthop isect];
                       included(isimp)=1;
                   end
               end
           end
           iynexthop=setdiff(iynexthop(:),iyneigh);
           iyneigh=[iyneigh; iynexthop];
           ihop=[ihop; repmat(Nhop,size(iynexthop))];
           Nleft=sum(~included);
           if(Nleft&&isempty(iynexthop))
               Nhop=Nhop+1;
               iynexthop=setdiff(yada(:),iyneigh);
               iyneigh =[iyneigh; iynexthop];
               ihop=[ihop; repmat(Nhop,size(iynexthop))];
               Nleft=0;
           else
               iynexthop=[];
           end
       end
       pwem(iy).iyneigh=iyneigh;
       pwem(iy).ihop=ihop;
       pwem(iy).Nneigh=length(pwem(iy).iyneigh);

       %iyneigh contains in order, this sample's id (0 hops away), the id 
       %of samples 1 hop away, the id of samples 2 hops away, etc. for as
       %many hops used in the definition of the neighborhood of the
       %mini-emulator
       
       if(1)
           %"new and improved" logic, in the process of being tested...
           %lmin=1/4 the radius of a hypersphere with the same "content"
           %(hypervolume) as the sum of the 1 hop simplices.  But to make 
           %this meaningful the different macro inputs must have comparable
           %scales, so multiply by the scale factors before computing the
           %content of the simplices and divide by the scale factors after
           %you've found the radius.  The radius of the 1 hop hypersphere
           %equals the average distance to adjacent points/simulations, we
           %want a "completely unsure" emulator to say it is most unsure
           %about the point 1/2 this distance away so say 1/2 radius = 2
           %correlation lengths (standard deviations) 
           %=> minium correlation length, lmin, = radius/4
           iy1hop=pwem(iy).iyneigh(find(pwem(iy).ihop<=1));
           content=0;
           for is=1:pwem(iy).Nsimpneigh
               iysimp=tess(pwem(iy).isimpneigh(is),:);
               if(isempty(setdiff(iysimp,iy1hop)))
                   content=content+simplex_content(y(iysimp,:).*SFM);
               end
           end
           pwem(iy).lmin=0.25*...
               (content/Sn*NdiminmacroY)^(1/NdiminmacroY)./scalefactor;       
       else 
           %"old" logic... 
           %but it hasn't yet been tried for the hazmap emulator so there
           %are no guarantees that this works well either
           iy1hop=pwem(iy).iyneigh(find(pwem(iy).ihop==1));    
           N1hop=length(iy1hop);
           pwem(iy).lmin=0.25*sqrt(...
               mean(((y(iy1hop,:)-repmat(y(iy1hop,:),N1hop,1)).*...
                    repmat(scalefactor,N1hop,1)).^2)...
               *NdiminmacroY)./scalefactor;       
       end
    end
    %save DEBUGME3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ifplot)       
       color=repmat('rmbcgykw',1,4);
       iy=samplenumber;
       ipts=tess(pwem(iy).isimpneigh,[1:3 1])';
       y1=reshape(y(ipts,1),4,pwem(iy).Nsimpneigh);
       y2=reshape(y(ipts,2),4,pwem(iy).Nsimpneigh);       
       figure; hold on;
       Nhop=max(pwem(iy).ihop);
       legendstr='legend(';
       for ihop=0:Nhop
           legendstr=sprintf('%s''%g hops'',',legendstr,ihop);
           jy=pwem(iy).iyneigh(find(ihop==pwem(iy).ihop));
           plot(y(jy,1),y(jy,2),['o' color(ihop+1)],'markerfacecolor',color(ihop+1),'markersize',7);
       end
       plot(y1,y2,'-k',y(iy,1),y(iy,2));
       legendstr=[legendstr '''location'',''best'');'];       
       hold off
       save DEBUGME5;
       title(sprintf('neighborhood of mini-emulator %g',samplenumber),'fontsize',14);
       axis normal;
       eval(legendstr);
    end
    
    fid=fopen(sprintf('build_mini_pwem_meta.%06g',samplenumber),'w');
    fprintf(fid,'additional file format lines=%g\n',4);

    fprintf(fid,'%%Ndiminmacro=4: #of macro-emulator input dimensions (log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
    fprintf(fid,'%%ymacro lmin: this mini-emulator''s minimum candidate correlation lengths for the macro variables: [Ndiminmacro] doubles\n');
    fprintf(fid,'%%Nneighmacro: #of simulations in the macro-emulator neighborhood of this simulation/mini-emulator: [1] integer\n');
    fprintf(fid,'%%{{(ihop): number of hops from macro center of mini-emulator: [1] integer},{iymacro: indices of simulations in the neighborhood of this mini-emulator, the first index in the list is the index/id of this simulation/mini-emulator, it''s also the filename''s suffix: [Nneighmacro] integers}}\n');
    
    fprintf(fid,sprintf('%%g\n%%g%s\n%g\n',...
        repmat(' %g',[1 NdiminmacroY-1])),NdiminmacroY,...
        pwem(samplenumber).lmin,pwem(samplenumber).Nneigh);
    fprintf(fid,'(%g) %g\n',...
        [pwem(samplenumber).ihop pwem(samplenumber).iyneigh]');
    fclose(fid);
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V=simplex_content(y)
    %the "content" of a simplex is its hypervolume, V.  See
    %http://mathworld.wolfram.com/Cayley-MengerDeterminant.html
    %y
    Ndim=size(y,2);
    Nnodes=size(y,1);
    if(Nnodes~=Ndim+1)
        error('Xminimetasimpcontbad','#of nodes in "simplex" (%g) ~= #of dimensions (%g) +1',Nnodes,Ndim);
    end
    
    B=double(~eye(Ndim+2)); %need the double(...) because of weird glitch 
    %in matlab's implicit typecasting which won't let B be treated as 
    %anything but a "logical" variable
    for i=1:Ndim
        for j=i+1:Nnodes
            yada=norm(y(i,:)-y(j,:));
            B(i+1,j+1)=yada;
            B(j+1,i+1)=yada;
            %disp(sprintf('i=%g j=%g B(i+1,j+1)=%g\n',i,j,yada));
        end
    end
    %B
    V=sqrt((-1)^(Ndim+1)/(2^Ndim*factorial(Ndim)^2)*det(B));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sn=hypersphere_Sn(Ndim)
    %the "content" (hypervolume), V, of a hypersphere is given by
    %V=(Sn*R^n)/n where n=Ndim, this function returns Sn. See
    %http://mathworld.wolfram.com/Hypersphere.html
    Sn=2*pi^(Ndim/2)/gamma(Ndim/2);
return
