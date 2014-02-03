function identify_mini_emulators_to_evaluate()
    ifplot=0; %1;
    tic;
    
    %disp('yada1');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read the resample inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid=fopen('macro_resamples.tmp','r');
    Nskip=sscanf(fgets(fid),'additional file format lines=%g',1);
    for i=1:Nskip
        fgets(fid);
    end
    NdiminmacroX=sscanf(fgets(fid),'%g',1);
    Nxmacro=sscanf(fgets(fid),'%g',1);
    xw=fscanf(fid,'%g',[NdiminmacroX+1 Nxmacro])';
    fclose(fid);

    %disp('yada2');
    
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

    if(NdiminmacroY~=NdiminmacroX)
        %need to display an error message
        return;
    end
        
    pwem=repmat(struct('Nneigh',0,'iyneigh',[],'isimpneigh',[]),Nymacro,1);

    %Nymacro
    %NdiminmacroY
    
    y=zeros(Nymacro,NdiminmacroY);
    %disp('yada3')
    for iy=1:Nymacro
        thisline=fgets(fid);
        yada=sscanf(thisline,sprintf('(%d) %s',iy,...
            repmat(' %g',1,NdiminmacroY+1)),NdiminmacroY+1);
        y(iy,:)=yada(1:NdiminmacroY);
        pwem(iy).isimpneigh=sscanf(thisline,sprintf('(%d)%s%s',iy,...
            sprintf(' %g',yada),repmat(' %g',1,yada(NdiminmacroY+1))),...
            yada(NdiminmacroY+1));
        pwem(iy).Nsimpneigh=yada(NdiminmacroY+1);
    end
%disp('yada4')
    
    Nsimpmacro=sscanf(fgets(fid),'%g',1);
    yada=sprintf('(%%g)%s\n',repmat(' %g',1,NdiminmacroY+1));
    tess=fscanf(fid,yada,[NdiminmacroY+2 Nsimpmacro])';
    fclose(fid);
%disp('yada4.5')

%size(tess)
    if(sum(abs(diff(tess(:,1))-1))||(tess(1,1)~=1))    
        error('need an error message');
        %if(any(tess(:,1)~=(1:Nsimpmacro)'))
        %need to display an error message
        return;
    end
    
    %disp('yada4.7');
    tess=tess(:,2:NdiminmacroY+2);
    
    %disp('yada5')
    if(ifplot)
        trimesh(tess,y(:,1),y(:,2));
        title('read in macro emulator','fontsize',14);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %this information is useful for building mini-emulators, not so much 
    %for identifying mini-emulators to evaluate
    %for iy=1:Nymacro
if(0)
    for iy=[]
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
       %mini-emulator, this is potentially useful information so since it 
       %didn't cost any extra storage (actually it saves a little by not
       %requiring iy0 to be stored) I built it into the data format.
    end
end
    %save DEBUGME3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %save DEVELOPEME
     
    %disp(sprintf('started tsearchn at time %g',toc));
    %figure; title(sprintf('started tsearchn at time %g',toc));
    %drawnow;
    %tic
    %[isimp,barycoord]=tsearchn(y,tess,xw(:,1:NdiminmacroY));
    %MATLABtime=toc
    
    %tic
    [isimp,barycoord]=mytsearchn(y,tess,xw(:,1:NdiminmacroY));
    %mytime=toc
    
    %ixwisnan=find(isnan(isimpmy));
    %Numisnan=numel(ixwisnan)
    %if(Numisnan)
    %   [isimp,barycoord]=tsearchn(y,tess,xw(ixwisnan,1:NdiminmacroY));
    %end
    
    %bob
    %disp(sprintf('finished tsearchn at time %g',toc));
    [isimp,isort]=sort(isimp);
    %isimp2=isimp;
    barycoord=barycoord(isort,:);
    xw=xw(isort,:); %sorting by tess makes it quicker/simpler to evaluate
    %the mini/macro-emulator and build the hazard map, all points within
    %a particular macro-emulator simplex are taken care of at the same 
    %time
    inan=find(isnan(isimp),1):length(isimp);
    inotnan=1:inan(1)-1;
    isimp=isimp(inotnan);
    barycoord=barycoord(inotnan,:);
    
    fid=fopen('macro_resample_assemble.inputs','w');
    fprintf(fid,'additional file format lines=%g\n',5);
   
    fprintf(fid,'%%Ndiminmacro=4: of macro-emulator input dimensions (log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
    fprintf(fid,'%%Nxmacroinside: the number of macro emulator (log10(volume [m^3]),Estart [UTME],Nstart [UTMN],BedFrictAng [deg])input points to evaluate: [1] integer\n');    
    fprintf(fid,'%%{{(log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [Ndiminmacro] doubles},{isimp: index of  macro-emulator simplex: [1] integer},{simulation/mini-emulator indices/ids: [Ndiminmacro+1] integers},{barycentric coordinates: [Ndiminmacro+1] doubles},{w: relative weight for hazmap assembly, sum(w)~=1 is ok: [1] double}}\n');
    fprintf(fid,'%%Nxmacrooutside: the number of randomly generated input points that can''t evaluate because they lie outside the convex hull of simulations\n');
    fprintf(fid,'%%{{(log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [Ndiminmacro] doubles},{w: relative weight for hazmap assembly: [1] double}}\n');
    
    fprintf(fid,'%d\n%d\n',NdiminmacroY,inan(1)-1);
    yada=sprintf('%s%%.10g\n',repmat('%.10g ',1,NdiminmacroY+1+(NdiminmacroY+1)*2));
    %save DEBUGME4;

    fprintf(fid,yada,[xw(inotnan,1:NdiminmacroY) isimp tess(isimp,:) barycoord xw(inotnan,NdiminmacroY+1)]');
    fprintf(fid,'%d\n',length(inan));

    yada=sprintf('%s%%.10g\n',repmat('%.10g ',1,NdiminmacroY));
    fprintf(fid,yada,xw(inan,:)');
    fclose(fid);
return