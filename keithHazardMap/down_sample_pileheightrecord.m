function down_sample_pileheightrecord(samplenumber)
    ifplot=0; %1
    ifwrite=1;
    if(ischar(samplenumber))
      samplenumber=str2num(samplenumber);
    end
    
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
    yada=fscanf(fid,'%g',[Ndimin Ny])';
    yada=yada(samplenumber,:);
    fclose(fid);
    
    log10vol   =yada(1);
    Direction  =yada(2);
    BedFrictAng=yada(3);
    IntFrictAng=yada(4);
    
    %fid=fopen(...
    %    sprintf('VolENstartBedSample.%06g',...
    %    samplenumber),'r');
    %fgetl(fid); fgetl(fid);
    %yada=fscanf(fid,'log10(Vol [m^3])=%g\nUTM East  Center=%g\nUTM North Center=%g\nBed Frict [deg] =%g',4);
    %fclose(fid);

    %log10vol     =yada(1);
    %UTMEastStart =yada(2);
    %UTMNorthStart=yada(3);
    %BedFrictAng  =yada(4);
    
    volume=10^log10vol; 
    GEOFLOW_TINY=volume^(1/3)/10000;

    filename=sprintf('pileheightrecord.%06g',samplenumber);
    fid=fopen(filename,'r');
    if(fid==-1)
        %need to do some kind of error code thing here
        return; 
    end
    fclose(fid);
    XYH=read_in_zgrid(filename);
    
    X=XYH(:,:,1);
    Y=XYH(:,:,2);
    H=XYH(:,:,3);
    clear XYZ s w;
    [Nx,Ny]=size(X);

    dx=diff(X(1:2,1));
    dy=diff(Y(1,1:2));
    TOOCLOSE=2*(dx+dy);

    ifnoH=(H<GEOFLOW_TINY);
    H(find(ifnoH))=0;
    WEIGHT=~ifnoH;
    deaddist=(sum(WEIGHT(:))*dx*dy/pi)^0.5;
    isd=ceil(0.1*deaddist/(0.5*(dx+dy))); isd=5;


    [garb,iMax]=max(H(:));
    xMax=X(iMax);
    yMax=Y(iMax);
    minmaxh=[GEOFLOW_TINY H(iMax)];
    CLev0=minmaxh(1); %contour level zero (lowest contour)
    CLev=CLev0*10.0.^((0:0.05:0.95)*log10(minmaxh(2)/CLev0));
    NCLev=length(CLev);

    H2=[zeros(Nx+2,1) [zeros(1,Ny); H; zeros(1,Ny)] zeros(Nx+2,1)];
    X2=repmat([X(1,1)-dx; X(:,1); X(Nx,1)+dx],[1    Ny+2]);
    Y2=repmat([Y(1,1)-dy  Y(1,:)  Y(1,Ny)+dy],[Nx+2 1   ]);

    if(ifplot)    
        [cc,hc]=contour(X2,Y2,H2,CLev,'-g');
        hold on;
    else
        %sx2=size(X2)
        %sy2=size(Y2)
        %sh2=size(H2)
        cc=contourc(X2(:,1)',Y2(1,:)',H2',CLev);
    end
    clear X2 Y2 H2;

    ixbase=[reshape(repmat([1 ceil(Nx/2) Nx]',[1 3]),[9 1]); ...
        reshape(repmat([ceil(Nx/4) ceil(0.75*Nx)]',[1 2]),[4 1])];
    iybase=[reshape(repmat([1 ceil(Ny/2) Ny] ,[3 1]),[9 1]); ...
        reshape(repmat([ceil(Ny/4) ceil(0.75*Ny)] ,[2 1]),[4 1])];
    ibase=unique((iybase-1)*Nx+ixbase);
    xbase=X(ibase);
    ybase=Y(ibase);

    itooclose=find(((xbase-xMax).^2+(ybase-yMax).^2).^0.5<=TOOCLOSE);
    ibase=ibase(setdiff(1:length(ibase),itooclose));    
    isofar=[iMax; ibase];
    Nsofar=length(isofar);
    NsofarOrig=Nsofar;
    xsofar=X(isofar);
    ysofar=Y(isofar);
    if(ifplot)
        plot(X(iMax),Y(iMax),'m*',X(ibase),Y(ibase),'mo');
    end
        
    ic=1;
    sizcc=size(cc);
    HULL=[];
    while(ic<sizcc(2))
        nc=cc(2,ic);  
        c=cc(:,ic+(1:nc));
        iCLev=find(cc(1,ic)==CLev);
        tooclose=((iCLev-1)/(NCLev-1)*(5-2)+2)*(dx+dy);
        TOOSMALL=tooclose^2;

        ich=convhull(c(1,:),c(2,:));
        HULL=[HULL c(:,ich)];
        [xx,yy,Numcut]=downsampleboundary(c(1,:),c(2,:),1.05);
        i=1:length(xx)-1;
        xmid=0.5*(xx(i)+xx(i+1));
        ymid=0.5*(yy(i)+yy(i+1));  
        parea=polyarea(xx,yy);
        Numleft=length(xx);
        cB=[xx; yy];  
  
        if(parea<=TOOSMALL)
            cB=mean(cB,2);
        else
            ds=diff(cB,[],2);
            magds=sum(ds.^2,1).^0.5;
            ikeep=setdiff(1:Numleft,setdiff(find(magds<TOOCLOSE),[1 Numleft]));
            cB=cB(:,ikeep);    
        end
        ixtmp=max(min(round((cB(1,:)-X(1,1))/dx+0.5),Nx),1);
        iytmp=max(min(round((cB(2,:)-Y(1,1))/dy+0.5),Ny),1);
  
        if(iCLev==1)
            Ntmp2=length(ixtmp);
            ixtmp2=zeros(size(ixtmp));
            iytmp2=ixtmp2;
            for itmp2=1:Ntmp2
                [ix,iy,ixmax,iymax]=find_nearest_boundary_and_max(X,Y,H,...
                ixtmp(itmp2),iytmp(itmp2),isd);
                ixtmp(itmp2)=ix;
                iytmp(itmp2)=iy;
                ixtmp2(itmp2)=ixmax;
                iytmp2(itmp2)=iymax;
            end
            ixtmp=[ixtmp2 ixtmp];
            iytmp=[iytmp2 iytmp];
        end
    
        itmp=unique((iytmp-1)*Nx+ixtmp);  
        xtmp=X(itmp);
        ytmp=Y(itmp);

        Ntmp=length(itmp);

        jtmp=[];
        for i=1:Ntmp  
            dist2=(xsofar-xtmp(i)).^2+(ysofar-ytmp(i)).^2;
            if(all(dist2>TOOSMALL))
                Nsofar=Nsofar+1;
                jtmp=[jtmp; itmp(i)];
                isofar(Nsofar)=itmp(i);
                xsofar(Nsofar)=xtmp(i);
                ysofar(Nsofar)=ytmp(i);
            end
        end
  
        if(ifplot)
            plot(X(jtmp),Y(jtmp),'r+');
        end
        ic=ic+nc+1;
    end
    TOOSMALL=TOOCLOSE^2;

    HULL=HULL(:,convhull(HULL(1,:),HULL(2,:)));
    xh=HULL(1,:);
    yh=HULL(2,:);

    disp('looking for WEIGHT');
    ifnoH=(H<GEOFLOW_TINY);
    WEIGHT=~ifnoH;
    deaddist=0.5*(sum(WEIGHT(:))*dx*dy/pi)^0.5;
    isd=ceil(deaddist/(0.5*(dx+dy)));
    [ix_yesh,iy_yesh]=find(WEIGHT);
%iyesh=(iy_yesh-1)*Nx+ix_yesh;
    ixrectminmax=max(min(...
        round([min(ix_yesh)-deaddist/dx max(ix_yesh)+deaddist/dx]),Nx),1);
    iyrectminmax=max(min(...
        round([min(iy_yesh)-deaddist/dy max(iy_yesh)+deaddist/dy]),Ny),1);
    WEIGHT=max(WEIGHT,0.05);


    for iy=iyrectminmax(1):iyrectminmax(2)  
        iiy=max(iy-isd,1):min(iy+isd,Ny);
        for ix=ixrectminmax(1):ixrectminmax(2)
            iix=max(ix-isd,1):min(ix+isd,Nx);
            hhh=H(iix,iiy);
            xxx=X(iix,iiy);
            yyy=Y(iix,iiy);
            if(H(ix,iy)<=0)
                iih=find(hhh);
                if(~isempty(iih))
                    dist=min((xxx(iih)-X(ix,iy)).^2+(yyy(iih)-Y(ix,iy)).^2)^0.5;
                    if((dist>0)&&(dist<=2*(dx+dy)))
                        WEIGHT(ix,iy)=2;
                    else
                        WEIGHT(ix,iy)=0.05+(1-dist/deaddist)*0.95;
                    end
                end
            end
        end
    end
    %disp('found WEIGHT');



    
    Nwant=max(Nsofar,1000-Nsofar);
    Nmore=0;
    while(Nmore<Nwant)
        %TOOSMALL
        trit=delaunay(xsofar,ysofar,{'QJ'})';
        %parea=polyarea(xsofar(trit),ysofar(trit));
        edgelen=(diff(xsofar(trit([1:3 1],:))).^2+...
                 diff(ysofar(trit([1:3 1],:))).^2).^0.5;
             pseudoparea=(3^0.5/4)*mean(edgelen,1).^2;  %the area of an 
             %equalateral triangle whose edge length is the mean edge 
             %length.
             [garb,imaxedg]=max(edgelen);
                   
             ipts=trit([imaxedg; mod(imaxedg,3)+1]+repmat(3*(0:size(trit,2)-1),[2 1]));
             xedg=mean(xsofar(ipts),1);
             yedg=mean(ysofar(ipts),1);
             ixedg=max(min(round((xedg-X(1,1))/dx+0.5),Nx),1);
             iyedg=max(min(round((yedg-Y(1,1))/dy+0.5),Ny),1);
             iedg=(iyedg-1)*Nx+ixedg;
             xedg=X(iedg);
             yedg=Y(iedg);
             %hedg=H(iedg);
             pseudoparea2=WEIGHT(iedg).*pseudoparea;
             [garb,isort]=sort(pseudoparea2,'descend');

             
             itmp=iedg(isort);
             xtmp=xedg(isort);
             ytmp=yedg(isort);
  
             Ntmp=length(itmp);
             dist2=min(((repmat(xsofar,[1 Ntmp])-repmat(xtmp,[Nsofar 1])).^2+...
                 (repmat(ysofar,[1 Ntmp])-repmat(ytmp,[Nsofar 1])).^2),[],1);
             iok=find(dist2>TOOSMALL); 
             Nok=min(numel(iok),5);
             itmp=unique(itmp(iok(1:Nok))');  
             Ntmp=length(itmp);  
             Nmore =Nmore +Ntmp;
             Nsofar=Nsofar+Ntmp;
             if(Ntmp==0)
                 Nwant=Nmore;
             end
             
             xtmp=X(itmp);
             ytmp=Y(itmp);
             if(ifplot)
                 plot(xtmp,ytmp,'m+');  
                 axis image;
                 drawnow;
             end

             isofar=[isofar; itmp];
             xsofar=[xsofar; xtmp];
             ysofar=[ysofar; ytmp];
    end

    x=X(isofar);
    y=Y(isofar);
    h=H(isofar);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the down sampled data file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %tess=delaunay(x,y,{'QJ'});
    %ikeeptess=find(any(h(tess)>=GEOFLOW_TINY,2));
    %tess=tess(ikeeptess,:);
    %ptstokeep=unique(tess(:));
    %inew=(1:numel(ptstokeep))';
    %inewfromiold=zeros(size(x));
    %inewfromiold(ptstokeep)=inew;
    %xkeep=x(ptstokeep); ykeep=y(ptstokeep); hkeep=h(ptstokeep);
    %tesskeep=inewfromiold(tess);
   
    
    
    if(ifwrite)        
        %save DEVELOPME;
        write_down_sampled_data(samplenumber,log10vol,Direction,BedFrictAng,IntFrictAng,x,y,h,minmaxh);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    if(~ifplot)
        return;
    end    
    H2=griddata(x,y,h,X,Y,'linear',{'QJ','QbB'});
    H2(find(isnan(H2)))=0;
    
    %rmslog=10^(mean(mean(log10(abs(H-H2)).^2))^0.5) %inf because log(0)=-inf
    %mean(mean((H-H2).^2))^0.5
    rms=mean(mean((H-H2).^2))^0.5


    %%%%%%%%%%%%%%%%%%%
    %figures
    %%%%%%%%%%%%%%%%%%%

    Nquads=(Nx-1)*(Ny-1);
    iX=shiftdim([1:Nx-1]'*ones(1,Ny-1),-1);
    iY=shiftdim(ones(Nx-1,1)*[1:Ny-1],-1);
    quads=zeros(4,Nx-1,Ny-1);
    quads(1,:,:)=(iY-1)*Nx+iX;
    quads(2,:,:)=(iY-1)*Nx+iX+1;
    quads(3,:,:)=(iY-0)*Nx+iX+1;
    quads(4,:,:)=(iY-0)*Nx+iX;
    quads=reshape(quads,[4 Nquads]);
    
    cmap=colormap;
    cmaplength=size(cmap,1);

    figure;
    iHzero=find(ifnoH);
    icolor=round((log10(H)-log10(minmaxh(1)))/diff(log10(minmaxh))*cmaplength);
    icolor(find(icolor>cmaplength))=cmaplength;
    icolor(find(icolor<1))=1;
    cvert=cmap(icolor,:);
    cvert(iHzero,:)=1;
    
    patch(X(quads),Y(quads),H(quads),reshape(cvert(quads,:),4,Nquads,3));
    shading interp;
    axis image;
    hcl=camlight;
    caxis(minmaxh);
    colorbar('vert','yscale','log');

    figure;
    iH2zero=find(~H2);
    icolor=round((log10(H2)-log10(minmaxh(1)))/diff(log10(minmaxh))*cmaplength);
    icolor(find(icolor>cmaplength))=cmaplength;
    icolor(find(icolor<1))=1;
    cvert2=cmap(icolor,:);
    cvert2(iH2zero,:)=1;

    patch(X(quads),Y(quads),H2(quads),reshape(cvert2(quads,:),4,Nquads,3));
    shading interp;
    axis image;
    hcl=camlight;
    caxis(minmaxh);
    colorbar('vert','yscale','log');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write the DOWN SAMPLED DATA FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_down_sampled_data(samplenumber,log10vol,Direction,BedFrictAng,IntFrictAng,X,Y,H,minmaxh)
%function write_down_sampled_data(samplenumber,log10vol,alpha,ysy,tess,minmaxh)
    iftoc=0;
    
    %save DEVELOPME;
    
    Nx        =size(X,1);    
    
    inzh      =find(H>=minmaxh(1));
    Nnzh      =numel(inzh);
    tess      =delaunay(X,Y,{'QJ'});
    Ntess     =size(tess,1);
    tess      =sortrows(sort(tess,2));
    itessneed =find(any(H(tess)>=minmaxh(1),2));
    Ntessneed =numel(itessneed);
    itessextra=setdiff((1:Ntess)',itessneed);
    tess      =tess([itessneed; itessextra],:);
    iptsneed  =unique(tess(1:Ntessneed,:));
    Nptsneed  =numel(iptsneed);
    
    for iter=1:200
        inewfromiold=repmat(-3,Nx,1);
        inewfromiold(iptsneed)=-2;
        inewfromiold(inzh)=-1;
        iptsneed=(1:Nptsneed)';
        inzh=(1:Nnzh)';
        
        temptess=reshape(tess',Ntess*3,1);
        inew1=0;
        inew2=Nnzh;
        inew3=Nptsneed;
        
        for i=1:Ntess*3
            j=temptess(i);
            %save DEBUGMEWTF_downsample1;
            switch(inewfromiold(j))
                case -1,
                    inew1=inew1+1;
                    inewfromiold(j)=inew1;
                case -2,
                    inew2=inew2+1;
                    inewfromiold(j)=inew2;
                case -3,
                    inew3=inew3+1;
                    inewfromiold(j)=inew3;
                otherwise,
            end
        end
        if(~((inew1==Nnzh)&&(inew2==Nptsneed)&&(inew3==Nx)))
            save DEBUGMEWTF_downsample2;
            bob;
        end
        X(inewfromiold)=X;
        Y(inewfromiold)=Y;
        H(inewfromiold)=H;
        tess =sortrows(sort(inewfromiold(tess),2));
    end
                        
    tessN=tess(1:Ntessneed,:);        
    ysy=[X(iptsneed) Y(iptsneed) H(iptsneed)];
    Ny=Nptsneed;    
    %I did something funky here, I label only those point that are the
    %nodes of simplices in which at least one of the nodes has non-zero
    %pileheight as being "needed" and the rest as being "extra" I had tried
    %to keep only the "needed" points but that really screwed up tsearch
    %and tsearchn, which needs the triangulation to be a delaunay
    %triangulation i.e. a convex hull, so I also need the "extra" points
    %but only for searching, I will only build and evaluate the "needed"
    %micro emulators
    

    if(0)
        %X=ysy(:,1); Y=ysy(:,2); H=ysy(:,3);
        cmap=colormap;
        cmaplength=size(cmap,1);
        
        iHzero=find(H<=minmaxh(1));
        icolor=min(max(real(round((log10(H)-log10(minmaxh(1)))/diff(log10(minmaxh))*cmaplength)),1),cmaplength);
        cvert=cmap(icolor,:);
        cvert(iHzero,:)=1;
        trit=tess';
        patch(X(trit),Y(trit),H(trit),reshape(cvert(trit,:),[size(trit) 3]));
        shading interp;
        axis image;
        hcl=camlight;
        caxis(minmaxh);
        colorbar('vert','yscale','log');
        save DEVELOPME;
    end
        
    %not really a whole piecewise emulator, just the neighborhood
    %definitions.
    pwem=repmat(...
        struct('Nneigh',0,'itemp1',[],'itemp2',[],'iyneigh',[],'ihop',[]),...
        Ny,1);
            
    for iy=1:Ny
        [i,j]   =find(tessN==iy);
        yadayada=unique(tessN(i,:)); 
        pwem(iy).itemp1=yadayada(:);
        pwem(iy).Nneigh=length(pwem(iy).itemp1);
        newbatchneigh=setdiff(pwem(iy).itemp1,iy);
        pwem(iy).iyneigh=[iy; newbatchneigh];
        pwem(iy).ihop=[0; ones(size(newbatchneigh))];
    end
    
    if(iftoc)
        disp(sprintf('found 1 hop spatial neighbors at t=%g [sec]',toc));
    end

    %>1 Nhops neighbors
    Nhop=2;
    for ihop=2:Nhop
        %ihop
        for iy=1:Ny
            for j=1:pwem(iy).Nneigh
                jy=pwem(iy).itemp1(j);
                pwem(iy).itemp2=[pwem(iy).itemp2; pwem(jy).iyneigh];
            end
            pwem(iy).itemp2=unique(pwem(iy).itemp2);
        end
        for iy=1:Ny
            newbatchneigh=setdiff(pwem(iy).itemp2,pwem(iy).iyneigh);
            pwem(iy).iyneigh=[pwem(iy).iyneigh; newbatchneigh];
            pwem(iy).ihop=[pwem(iy).ihop; repmat(ihop,size(newbatchneigh))];
            pwem(iy).itemp2=[];
        end
        if(iftoc)
            disp(sprintf('found %g hop spatial neighbors at t=%g [sec]',ihop,toc));
        end
    end
    %disp('neighborhoods have been built');
    for iy=1:Ny
        pwem(iy).itemp1=[];
        pwem(iy).Nneigh=length(pwem(iy).iyneigh);                
    end

    fid=fopen(sprintf('down_sampled_data.%06g',samplenumber),'w');
    fprintf(fid,'additional file format lines=%g\n',7);
    fprintf(fid,'%%Ndiminmacro=4: of macro-emulator input dimensions (log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
    fprintf(fid,'%%ymacro=(log10(volume [m^3]),Direction [deg CC from East],BedFrictAng [deg],IntFrictAng [deg]): values of "uncertain inputs" for the particular simulation: [Ndiminmacro] doubles\n');
    fprintf(fid,'%%Ny Nycover: #of point and #of points needed to cover the flow: [2] integers\n');
    fprintf(fid,'%%[Nycover] lines of point info: each contains {{(iy): pt id: [1] integer},{y=(east, north, h) [3] doubles},{Nneigh: # of pts in neighborhood: [1] integer},{iyneigh: indices of pts in neighborhood: [Nneigh] integers},{neighhops: #number of hops each pt in neighborhood is away from iy: [Nneigh] integers}}\n');
    fprintf(fid,'%%[Ny-Nycover] lines of limited point info: each contains {{(iy): pt id: [1] integer},{y=(east,north): we know h=0 for these: [2] doubles}}\n');
    
    fprintf(fid,'%%{Ntri Ntricover}: #of triangles in east-north tessellation and the #of triangles needed to cover the flow: [2] integers\n');
    fprintf(fid,'%%tri: triangle node indices (into y): [Ntri 3] array of doubles\n');
    
    %save DEBUGMEWTF_downsample3;
    fprintf(fid,'%g\n%.10g %.10g %.10g %.10g\n%d %d\n',4,log10vol,Direction,BedFrictAng,IntFrictAng,Nx,Ny);
    for iy=1:Ny
        yada=sprintf('(%g) %.2f %.2f %g %g%s\n',iy,ysy(iy,:),pwem(iy).Nneigh,repmat(' %g',1,2*pwem(iy).Nneigh));
        fprintf(fid,yada,pwem(iy).iyneigh,pwem(iy).ihop);
    end
    iy=(Nptsneed+1:Nx)';
    fprintf(fid,'(%g) %.2f %.2f\n',[iy X(iy) Y(iy)]');

    fprintf(fid,'%g %g\n',Ntess,Ntessneed);
    fprintf(fid,'%g %g %g\n',tess');
        
    fclose(fid);    
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xyz=read_in_zgrid(filename)

    fid=fopen(filename,'r');
    yada=fscanf(fid,'Nx=%g: X={%g,%g}\nNy=%g: Y={%g,%g',[3 2]);
    fgets(fid); fgets(fid); %'}\nElevation=\n' OR '}\nPileheight=\n'
    xyz=zeros(yada(1,1),yada(1,2),3);
    xyz(:,:,1)=((2*(0:(yada(1,1)-1))+0.5)/(2*yada(1,1))*(yada(3,1)-yada(2,1))+yada(2,1))'*ones(1,yada(1,2)); %';
    xyz(:,:,2)=ones(yada(1,1),1)*(2*(0:(yada(1,2)-1))+0.5)/(2*yada(1,2))*(yada(3,2)-yada(2,2))+yada(2,2);
    xyz(:,:,3)=fscanf(fid,'%g',[yada(1,1) yada(1,2)]);
    fclose(fid);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,Numcut]=downsampleboundary(x,y,crit)
    %crit=0.98 is recommended
    Numcut=0;

    sizx=size(x);
    [xdimmax,ixdimmax]=max(sizx);
    sizy=size(y);
    [ydimmax,iydimmax]=max(sizy);
    if((prod(sizx)~=xdimmax)||(prod(sizy)~=ydimmax))
        error('x and y are not 1 dimensional arrays')
    elseif((xdimmax~=ydimmax)||(ixdimmax~=iydimmax))
        error('x and y are of different sizes');
    end
    N=xdimmax;
    %if(N<=5)
    %    return;
    %end
    if(~((x(1)==x(N))&&(y(1)==y(N))))
        N=N+1;
        x(N)=x(1);
        y(N)=y(1);
    end
    
    xy=[x(:) y(:)]';
    icut=0;
    while(~isempty(icut))
        ds=diff(xy,[],2);
        magds=sum(ds.^2,1).^0.5;
        ds=ds./repmat(magds,[2 1]);
        i=[N-1 1:N-2];
        dsdotds=sum(ds.*ds(:,i),1);
        dstemp=dsdotds./(0.5*(magds+magds(i)))*mean(magds);
        ilarge=find(dstemp>crit);
        ilarge2=ilarge(1:2:length(ilarge));
        icut=setdiff(ilarge,ilarge2+1); %don't delete adjacent points in 
        %the same round
        i=setdiff(i,icut);    
        i(length(i)+1)=i(1);
        N=length(i);
        xy=xy(:,i);
        Numcut=Numcut+length(icut);
    end
    sizx(:)=1;
    sizx(ixdimmax)=N;
    x=reshape(xy(1,:),sizx);
    y=reshape(xy(2,:),sizx);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ix,iy,ixmax,iymax]=find_nearest_boundary_and_max(...
    X,Y,H,ix,iy,isd)     
    [Nx,Ny]=size(H);
    %X,Y,H are for large map
    %ix,iy identify a point for whom you want to find the nearest boundary
    %of the flow and a maximum nearby flow depth
    %isd tells you how many squares away from ix,iy you want to search
    
    iix=max(ix-isd,1):min(ix+isd,Nx);
    iiy=max(iy-isd,1):min(iy+isd,Ny);
    Iix=repmat(iix',[1 length(iiy)]);
    Iiy=repmat(iiy ,[length(iix) 1]);
    hhh=H(iix,iiy);
    xxx=X(iix,iiy);
    yyy=Y(iix,iiy); 
    if(H(ix,iy)<=0)
        iih=find(hhh);
        if(~isempty(iih))
            [garb,imin]=min((xxx(iih)-X(ix,iy)).^2+(yyy(iih)-Y(ix,iy)).^2);
            iih=iih(imin);
            ix=Iix(iih);
            iy=Iiy(iih);
        end
    end
    if(H(ix,iy)>0)
        iih=find(~hhh);
        if(~isempty(iih))
            [garb,imin]=min((xxx(iih)-X(ix,iy)).^2+(yyy(iih)-Y(ix,iy)).^2);
            iih=iih(imin);
            ix=Iix(iih);
            iy=Iiy(iih);
        end
        [garb,imax]=max(hhh(:));
        ixmax=Iix(imax);
        iymax=Iiy(imax);
    else
        ixmax=ix;
        iymax=iy;
    end
return