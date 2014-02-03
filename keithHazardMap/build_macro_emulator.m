function build_macro_emulator()
    ifplot=0; %1;
    fid=fopen('uncertain_input_list.txt','r');
    if(fid==-1)
        %send an error code or message
        return;
    end
    Nskip=sscanf(fgets(fid),'additional file format lines=%g',1);
    for i=1:Nskip;
        fgets(fid);
    end
    Ndimin=fscanf(fid,'%g',1);
    Ny=fscanf(fid,'%g',1);
    y=fscanf(fid,'%g',[Ndimin Ny])';
    fclose(fid);
  
    fid=fopen('macro_emulator.pwem','w');
    fprintf(fid,'additional file format lines=%g\n',5);

    fprintf(fid,'%%Ndiminmacro=4: of macro-emulator input dimensions (log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
    fprintf(fid,'%%Nymacro: #number of simulations/mini-emulators: [1] integer\n');
    fprintf(fid,'%%[Nymacro] lines of simulation/mini-emulator info: each contains {{(iy): simulation id: [1] integer},{y=(log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg]): [Ndiminmacro] doubles},{Nsimpneigh: # of simplices in neighborhood: [1] integer},{isimpneigh: indices of simplices in neighborhood: [Nsimpneigh] integers}}\n');
    fprintf(fid,'%%Nsimpmacro: #number of simplices (in this case triangles) in macro-emulator: [1] integer\n');
    fprintf(fid,'%%tess:  (index of the simplex) followed by indices of simulation/min-emulator nodes in each simplex: a [Nsimpmacro Ndiminmacro+1] array of integers\n');
  
    fprintf(fid,'%d\n%d\n',Ndimin,Ny);
    
    tess=sortrows(sort(delaunayn(y,{'QJ','QbB'}),2)); %sorting within rows
    %then row sorting puting lowest node numbers first (lowest node number
    %in simplex counts the most, followed by the second lowest, etc.) means
    %that running simulations in the order of most importat ones first,
    %places the most important simplexes first in the list.  If adaptivity
    %is introduced, then adding points in important regions after the first
    %batch will not "hurt" the priority of the most important simplices
    %much because the lowest node numbers count most in the sorting scheme.
    %an added benefit is it puts the simplices in an order which "makes
    %sense" to human readers of the data file
  
    pwem=repmat(struct('Nneigh',0,'itemp1',[],'itemp2',[],...
      'iyneigh',[],'itemp3',[],'isimpneigh',[]),Ny,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %build neighborhoods for mini emulators 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %1 hop neighbors
  
    for iy=1:Ny
        [i,j]          =find(tess==iy);
        pwem(iy).isimpneigh=[pwem(iy).isimpneigh; i(:)];
        yadayada=unique(tess(i,:));
        pwem(iy).itemp1=yadayada(:);
        pwem(iy).Nneigh=length(pwem(iy).itemp1);
        newbatchneigh=setdiff(pwem(iy).itemp1,iy);
        pwem(iy).iyneigh=[iy; newbatchneigh];
    end 
  
    %>1 Nhops neighbors
    Nhop=1;
    for ihop=2:Nhop
        %ihop
        for iy=1:Ny
            for j=1:pwem(iy).Nneigh
                jy=pwem(iy).itemp1(j);
                pwem(iy).itemp2=[pwem(iy).itemp2; pwem(jy).iyneigh];
                pwem(iy).itemp3=[pwem(iy).itemp3; pwem(jy).isimpneigh];
            end
            pwem(iy).itemp2=unique(pwem(iy).itemp2); 
            pwem(iy).itemp3=unique(pwem(iy).itemp3);
        end
        for iy=1:Ny
            newbatchneigh=setdiff(pwem(iy).itemp2,pwem(iy).iyneigh);
            pwem(iy).iyneigh=[pwem(iy).iyneigh; newbatchneigh];
            pwem(iy).itemp2=[];
            newbatchneigh=setdiff(pwem(iy).itemp3,pwem(iy).isimpneigh);
            pwem(iy).isimpneigh=[pwem(iy).isimpneigh; newbatchneigh];
            pwem(iy).itemp3=[];
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iy=1:Ny
        yada=tess(pwem(iy).isimpneigh,:);
        if(~isempty(setxor(pwem(iy).iyneigh,yada(:))))%sanity check            
            iy
            iyneigh1=unique(pwem(iy).iyneigh)
            iyneigh2=unique(tess(pwem(iy).isimpneigh,:))
            return;
        end
        fprintf(fid,'(%d)%s %d%s\n',iy,sprintf(' %.10g',y(iy,:)),...
            length(pwem(iy).isimpneigh),sprintf(' %d',pwem(iy).isimpneigh));
    end
    
    Nsimp=size(tess,1);
    fprintf(fid,'%d\n',Nsimp);
    yada=sprintf('(%%d)%s\n',repmat(' %d',1,Ndimin+1));
    fprintf(fid,yada,[1:Nsimp; tess']);
    %for isimp=1:Nsimp        
    %    fprintf(fid,'(%d)%s\n',isimp,sprintf(' %d',tess(isimp,:)));
    %end
    fclose(fid);
    
    if(ifplot)
        trimesh(tess,y(:,1),y(:,2))
    end

return;
