function make_montserrat_take2_pwem_dem(volume,direction,BEDFRICTANG,INTFRICTANG,samplenumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%volume      = volume of dome collapse in [m^3]
%
%direction   = preferred direction measured in degrees counter clockwise
%              from east  
%samplenumber= an integer used to distinguisg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ifplot=0; %1;
    if(ischar(volume))
      volume=str2double(volume);
    end
    if(ischar(direction))
      direction=str2double(direction);
    end
    if(ischar(BEDFRICTANG))
      BEDFRICTANG=str2double(BEDFRICTANG);
    end
    if(ischar(INTFRICTANG))
      INTFRICTANG=str2double(INTFRICTANG);
    end    
    if(ischar(samplenumber))
      samplenumber=str2num(samplenumber);
    end

    disp(sprintf('volume=%g direction=%g bedfrictang=%g intfrictang=%g samplenumber=%g',volume,direction,BEDFRICTANG,INTFRICTANG,samplenumber));
       

    [errorcode,filenamespresent]=system('dir');
    if(isempty(findstr('QUICKSTART.mat',filenamespresent)));
        %only need to do this once to generate a .mat file then save some 
        %time by loading the .mat file
        fid=fopen('mosaic.asc','r');

        NX=fscanf(fid,'ncols %g',1);
        NY=fscanf(fid,'\nnrows %g',1);
        xllc=fscanf(fid,'\nxllcorner %g',1);
        yllc=fscanf(fid,'\nyllcorner %g',1);
        DX=fscanf(fid,'\ncellsize %g',1); DY=DX;
        ndv=fscanf(fid,'\nNODATA_value %g',1);
        z=fscanf(fid,'%g',[NX NY]);
        fclose(fid);

        xminmax=xllc+[0 DX*(NX-1)];
        yminmax=yllc+[0 DY*(NY-1)];
        
        x=(xminmax(1):DX:xminmax(2))'*ones(1,NY);
        y=ones(NX,1)*(yminmax(2):-DY:yminmax(1));

        i=find(z==ndv);
        ii=find(z~=ndv);
        minz=min(z(ii));
        z(i)=minz;
        z(find(isnan(z)))=minz;

        z=fixslope(x,y,z,60,10000);

        save QUICKSTART x y z NX NY DX DY xminmax yminmax;
    else
        load QUICKSTART;
    end
    
    gisname=sprintf('montserrat_take2_pwem_%06g',samplenumber);
    eval('!mkdir grass5');
    eval('!mkdir grass5/montserrat_take2_pwem');
    eval(['!mkdir grass5/montserrat_take2_pwem/' gisname]);
    eval(['!mkdir grass5/montserrat_take2_pwem/' gisname '/cellhd']);
    eval(['!mkdir grass5/montserrat_take2_pwem/' gisname '/fcell']);

    
    tic;
    %direction
    z=alter_dem_and_gen_pile(x,y,z,volume,direction,gisname,BEDFRICTANG,INTFRICTANG,ifplot);

    
    fid=fopen(['grass5/montserrat_take2_pwem/' gisname '/fcell/' gisname],'wb');
    fwrite(fid,z,'single',0,'ieee-be');
    fclose(fid);

    fid=fopen(['grass5/montserrat_take2_pwem/' gisname '/cellhd/' gisname],'w');
    fprintf(fid,'proj:       %g\n', 0);
    fprintf(fid,'zone:       %g\n', 0);
    fprintf(fid,'north:      %g\n', yminmax(2));
    fprintf(fid,'south:      %g\n', yminmax(1));
    fprintf(fid,'east:       %g\n', xminmax(2));
    fprintf(fid,'west:       %g\n', xminmax(1));
    fprintf(fid,'cols:       %g\n',NX);
    fprintf(fid,'rows:       %g\n',NY);
    fprintf(fid,'e-w resol:  %g\n',DX);
    fprintf(fid,'n-s resol:  %g\n',DY);
    fprintf(fid,'format:     %g\n',-1);%float =-1
    fprintf(fid,'compressed: %g\n', 0);%no=0
    fclose(fid);
    toc

return;

function z=alter_dem_and_gen_pile(x,y,z,...
    pilevolume,cutpointingdown,gisname,BEDFRICTANG,INTFRICTANG,ifplot)
%parameters you set start here, these values have been "approved of" by
%Prof Eliza Calder. these angles are in degrees
conesteepangle=45; %cone is steeper on one end, this is the angle the 
%steeper end of the cone makes with the horizontal (positive only)
coneshallowangle=33; %cone is less steep on the opposite side, this is the
%angle at the least steep side makes with horizontal (positive only)
steepsidepointing=180; %the steep side is this many degrees counter 
%clockwise from east;
%cutangle=30; %angle at which dome is "cut." Above the cut plane is pile,
%below the cut plane is DEM (positive only);
cutangle=15; %angle at which dome is "cut." Above the cut plane is pile,
%below the cut plane is DEM (positive only);
%cutpointingdown=dirvol(1); %the downslope direction of the cut plane is this many 
%degrees clockwise from east
%domemaxheight=1050; %meters above sea level
coneratdomemaxheight=350; %radius in meters of "circular" top of "uncut" 
%DEM dome... if the dome was entirely represented by a DEM cone, the top
%would be elliptical in shape, but have the same area as a circle with a
%radius of this many meters.
centerdometopeast = 381000; %DEM east  coordinate of uncut dome top center
centerdometopnorth=1847100; %DEM north coordinate of uncut dome top center
%pilevolume=dirvol(2);
domemaxheight=max(1050,250+(2/pi*pilevolume)^(1/3)); %meters above sea level


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%quantities I compute start here, you don't change anything below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deg2rad=pi/180;
    Ntheta=720*16;
    theta=...  %angular position of point on unrotated cone in radians
        (0:Ntheta)/Ntheta*2*pi;
    gamma=...  %cone spread angle in radians
        (180-(conesteepangle+coneshallowangle))/2*deg2rad;
    alpha1=... %tilt angle in radians
        (conesteepangle-coneshallowangle)/2*deg2rad;
    alpha2=...
        steepsidepointing*deg2rad;
    beta1=...
        cutangle*deg2rad;
    beta2=...
        cutpointingdown*deg2rad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find ztop xoffset yoffset zoffset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ztemp=-(0:512)/512*1.125*coneratdomemaxheight/sin(gamma);
xtemp1=... %end of downslope edge of uncut top
    ztemp*(-sin(gamma)*cos(alpha1)/(sin(gamma)*sin(alpha1)+cos(alpha1))...
    +sin(alpha1)/(sin(gamma)*sin(alpha1)+cos(alpha1)));

xtemp2=... %end of upslope edge of uncut top
    ztemp*(sin(gamma)*cos(alpha1)/(-sin(gamma)*sin(alpha1)+cos(alpha1))...
    +sin(alpha1)/(-sin(gamma)*sin(alpha1)+cos(alpha1)));

xtemp3=... %focus of ellipse  
    ztemp*sin(alpha1)/cos(alpha1);
    
rmajortop=abs(xtemp1-xtemp2)/2;
xcentop=(xtemp1+xtemp2)/2;
deltax=abs(xcentop-xtemp3);
deltay=-ztemp*sin(gamma)/cos(alpha1);
rminortop=(0.25*(((2*deltax).^2+deltay.^2).^0.5+deltay).^2-deltax.^2).^0.5;
[garb,ibest]=min((rmajortop.*rminortop-coneratdomemaxheight^2).^0.5);
ztop=ztemp(ibest);
zoffset=-ztop+domemaxheight;
xoffset=xcentop(ibest)*cos(alpha2)+centerdometopeast;
yoffset=xcentop(ibest)*sin(alpha2)+centerdometopnorth;
clear ztemp xtemp1 xtemp2 xtemp3 rmajortop xcentop garb ibest deltax...
    deltay rminortop;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find zcut, pile height, pile major and minor radii and the counter
%clockwise from east angle of the pile's major radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp1=...
    sin(beta1)*cos(alpha2-beta2)*sin(gamma)*cos(theta)*cos(alpha1)...
    -sin(beta1)*cos(alpha2-beta2)*sin(alpha1)...
    -sin(beta1)*sin(gamma)*sin(theta)*sin(alpha2-beta2)...
    -cos(beta1)*sin(gamma)*cos(theta)*sin(alpha1)-cos(beta1)*cos(alpha1);

z4=-1;
x4=z4*...
    (-cos(beta1)*cos(alpha2-beta2)*sin(alpha1) ...
    -cos(beta1)*sin(gamma)*sin(theta)*sin(alpha2-beta2)...
    +cos(beta1)*cos(alpha2-beta2)*sin(gamma)*cos(theta)*cos(alpha1)...
    +sin(beta1)*sin(gamma)*cos(theta)*sin(alpha1)+sin(beta1)*cos(alpha1))...
    ./temp1;
y4=z4*...
    (sin(alpha2-beta2)*sin(gamma)*cos(theta)*cos(alpha1)...
    -sin(alpha1)*sin(alpha2-beta2)...
    +sin(gamma)*sin(theta)*cos(alpha2-beta2))...
    ./temp1;

x4focus=z4*...
    (-cos(beta1)*cos(alpha2-beta2)*sin(alpha1)+sin(beta1)*cos(alpha1))...
    /(-sin(beta1)*cos(alpha2-beta2)*sin(alpha1)-cos(beta1)*cos(alpha1));
y4focus=z4*...
    -(sin(alpha1)*sin(alpha2-beta2))...
    /(-sin(beta1)*cos(alpha2-beta2)*sin(alpha1)-cos(beta1)*cos(alpha1));

[garb,ithetamajor1]=max((x4-x4focus).^2+(y4-y4focus).^2);

ithetamajor2=mod(ithetamajor1-1+Ntheta/2,Ntheta)+1;

xtemp1=x4(ithetamajor1);
ytemp1=y4(ithetamajor1);
xtemp2=x4(ithetamajor2);
ytemp2=y4(ithetamajor2);

%maxy4=max(y4)
%r4maj=0.5*((xtemp1-xtemp2)^2+(ytemp1-ytemp2)^2)^0.5;
x4cen=(xtemp1+xtemp2)/2;
y4cen=(ytemp1+ytemp2)/2;
temp=(x4-x4cen).^2+(y4-y4cen).^2;
[garb,ithetaminor1]=min(temp);
if(mod(round(theta(ithetaminor1)-theta(ithetamajor1))/deg2rad+360,360)>180) 
  temp(ithetaminor1)=inf;
  [garb,ithetaminor1]=min(temp);
end
%r4min=temp(ithetaminor1)^0.5;

%clear xtemp1 ytemp1 xtemp2 ytemp2 r4maj r4min;
clear xtemp1 ytemp1 xtemp2 ytemp2;

ithetamajmin=[ithetamajor1 ithetaminor1];

x4mmc=[x4(ithetamajmin) x4cen x4focus];
y4mmc=[y4(ithetamajmin) y4cen y4focus];

x5mmc=x4mmc*cos(-beta1)-z4*sin(-beta1);
z5mmc=x4mmc*sin(-beta1)+z4*cos(-beta1);
y5mmc=y4mmc;

x6mmc=x5mmc*cos(beta2)-y5mmc*sin(beta2);
y6mmc=x5mmc*sin(beta2)+y5mmc*cos(beta2);
z6mmc=z5mmc;

i=1:4;
dR6mm=[x6mmc(i)-x6mmc(3); y6mmc(i)-y6mmc(3); z6mmc(i)-z6mmc(3)];
R6mm=sum(dR6mm(1:2,:).^2,1).^0.5; %scale by -z4 and enter in titan as major minor radius

Nrnd=128;
rnd=(0:Nrnd)'/Nrnd; %r non dimensional 0 to 1
rndct=rnd*cos(theta);
rndst=rnd*sin(theta);

phdh=(1-rnd.^2)*ones(size(theta));
phdh=phdh(:);
z6dz4=z6mmc(3)+rndct*dR6mm(3,1)+rndst*dR6mm(3,2);
z6dz4=z6dz4(:);

phdh=1; %(1-(R6mm(4)/R6mm(1))^2);
z6dz4=z6mmc(3); %z6mmc(4);

%phdh=(1-(R6mm(4)/R6mm(1))^2);
%z6dz4=z6mmc(4);

%max(z6dz4*-z4cut+phdh*h)=domemaxheight-zoffset; %-z4cut>0
%pi/2*h*R6mm(1)*R6mm(2)*z4cut^2=volume;

const1=domemaxheight-zoffset;
%const2=(2/(pi*R6mm(1)*R6mm(2))*pilevolume).^0.5; %for height newton
const2=(2/(pi*R6mm(1)*R6mm(2))*pilevolume); %for absz4cut newton

%absz4cut=abs(phdh*const2/z6dz4)^(1/3) %*sign(phdh*const2/z6dz4)

if(1)
    %h=(pilevolume*2/pi)^(1/3)*ones(size(phdh));
    h=(pilevolume*2/pi)^(1/3);
    absz4cut=(const2/h)^0.5;    
  


    for inewt=1:15
        %f=const2*z6dz4.*h.^-0.5+phdh.*h-const1;
        %df=-0.5*const2*z6dz4.*h.^-1.5+phdh;
        %d2f=0.75*const2*z6dz4.*h.^-2.5;
        f=z6dz4*absz4cut+phdh*const2/absz4cut^2-const1;
        [f,imax]=max(f);
        df =z6dz4(imax)-2*phdh(imax)*const2/absz4cut^3;
        %d2f=6*phdh(imax)*const2/absz4cut^4;
        %g=f.^2;
        %dg=2*f.*df;
        %d2g=2*(f.*d2f+df.^2);
        %dh=-f./df;
        %dh=-dg./abs(d2g);
        %dabsz4cut=-dg/abs(d2g);
        dabsz4cut=-f/df;
        %h=max(h+dh,0);
        absz4cut=max(absz4cut+dabsz4cut,1);
    end
end

if(pilevolume>0)
    pileheight=const2/absz4cut^2;
else
    pileheight=0;
end

pilemajorrad=R6mm(1)*absz4cut;
pileminorrad=R6mm(2)*absz4cut;

%pilevolume
pilevolume2=pi/2*pileheight*pilemajorrad*pileminorrad;
%cutangle
%cutpointingdown

xpilecenter=x6mmc(3)*absz4cut+xoffset;
ypilecenter=y6mmc(3)*absz4cut+yoffset;
zpilecenter=z6mmc(3)*absz4cut+zoffset;
pilemajoraxispointing=atan2(dR6mm(2,1),dR6mm(1,1))/deg2rad;

topofcone=max(absz4cut*z6dz4+pileheight*phdh)+zoffset;

%gisname1=sprintf('dirvol.%06d',isample);

fid=fopen(['grass5/montserrat_take2_pwem/' gisname '/' gisname '.pileprops'],'w');
%fid
fprintf(fid,'log10(volume [m^3])                        =%20.14g\n',log10(pilevolume2));
fprintf(fid,'cutangle                                   =%20.14g [deg]\n',cutangle);
fprintf(fid,'cut pointing down direction [cc from east] =%20.14g [deg]\n',cutpointingdown);
fprintf(fid,'hpilecenter                                =%20.14g [m]\n',pileheight);
fprintf(fid,'xpilecenter                                =%20.14g [East]\n',xpilecenter);
fprintf(fid,'ypilecenter                                =%20.14g [North]\n',ypilecenter);
fprintf(fid,'zpilecenter                                =%20.14g [m]\n',zpilecenter);
fprintf(fid,'major radius                               =%20.14g [m]\n',pilemajorrad);
fprintf(fid,'minor radius                               =%20.14g [m]\n',pileminorrad);
fprintf(fid,'major axis direction [cc from east]        =%20.14g [deg]\n',pilemajoraxispointing);
fprintf(fid,'bed friction angle                         =%20.14g [deg]\n',BEDFRICTANG);
fprintf(fid,'internal friction angle                    =%20.14g [deg]\n',INTFRICTANG);
fclose(fid);

%bob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add cone to DEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ntheta=360;
theta=...  %angular position of point on unrotated cone in radians
    (0:Ntheta)/Ntheta*2*pi;


Z1=-(0:10:2000+domemaxheight)';
X1=-Z1*sin(gamma)*cos(theta);
Y1=-Z1*sin(gamma)*sin(theta);
Z1=Z1*ones(size(theta));
z4=-absz4cut;
z1=...
    -z4./(sin(beta1)*cos(alpha2-beta2)*sin(gamma)*cos(theta)*cos(alpha1)...
    -sin(beta1)*cos(alpha2-beta2)*sin(alpha1)...
    -sin(beta1)*sin(gamma)*sin(theta)*sin(alpha2-beta2)...
    -cos(beta1)*sin(gamma)*cos(theta)*sin(alpha1)-cos(beta1)*cos(alpha1));
x1=-z1*sin(gamma).*cos(theta);
y1=-z1*sin(gamma).*sin(theta);


X2=X1*cos(-alpha1)-Z1*sin(-alpha1); 
Y2=Y1; 
Z2=min(X1*sin(-alpha1)+Z1*cos(-alpha1),ztop);
x2=x1*cos(-alpha1)-z1*sin(-alpha1);
y2=y1;
z2=x1*sin(-alpha1)+z1*cos(-alpha1);
clear X1 Y1 Z1 x1 y1 z1;

X3=X2*cos(alpha2-beta2)-Y2*sin(alpha2-beta2); 
Y3=X2*sin(alpha2-beta2)+Y2*cos(alpha2-beta2); 
Z3=Z2; 
x3=x2*cos(alpha2-beta2)-y2*sin(alpha2-beta2);
y3=x2*sin(alpha2-beta2)+y2*cos(alpha2-beta2);
z3=z2;
clear X2 Y2 Z2 x2 y2 z2;

X4=X3*cos(beta1)-Z3*sin(beta1); 
Y4=Y3;
Z4=min(X3*sin(beta1)+Z3*cos(beta1),-absz4cut);
x4=x3*cos(beta1)-z3*sin(beta1); 
y4=y3;
%z4=x3*sin(beta1)+z3*cos(beta1)=-absz4cut;
z4=-absz4cut*ones(size(z3));
clear X3 Y3 Z3 x3 y3 z3;

X5=X4*cos(-beta1)-Z4*sin(-beta1);
Y5=Y4; 
Z5=X4*sin(-beta1)+Z4*cos(-beta1);
x5=x4*cos(-beta1)-z4*sin(-beta1);
y5=y4; 
z5=x4*sin(-beta1)+z4*cos(-beta1);
clear X4 Y4 Z4 x4 y4 z4;

X6=X5*cos(beta2)-Y5*sin(beta2)+xoffset;
Y6=X5*sin(beta2)+Y5*cos(beta2)+yoffset;
Z6=Z5+zoffset;
x6=x5*cos(beta2)-y5*sin(beta2)+xoffset;
y6=x5*sin(beta2)+y5*cos(beta2)+yoffset;
z6=z5+zoffset;
clear X5 Y5 Z5 x5 y5 z5;

%surf(X6,Y6,Z6); shading flat; axis image; colormap white; clh=camlight; hold on
%plot3(x6,y6,z6,'-r');
%view(0,90);
%bob;

xx6=x6-xpilecenter; yy6=y6-ypilecenter;
rr6=(xx6.^2+yy6.^2).^0.5;
tt6=atan2(yy6,xx6);
%XX6=X6-xpilecenter; YY6=Y6-ypilecenter;
%RR6=(XX6.^2+YY6.^2).^0.5;
%TT6=atan2(XX6,YY6);

xx=x-xpilecenter; yy=y-ypilecenter;
rr=(xx.^2+yy.^2).^0.5;
tt=atan2(yy,xx);
%TT=tt-TT6(1);
tt=tt-tt6(1);

i=find(tt<0);    tt(i)=tt(i)+2*pi;
i=find(tt<0);    tt(i)=tt(i)+2*pi; 
i=find(tt>2*pi); tt(i)=tt(i)-2*pi;
i=find(tt>2*pi); tt(i)=tt(i)-2*pi;

%i=find(TT<0);    TT(i)=TT(i)+2*pi;
%i=find(TT<0);    TT(i)=TT(i)+2*pi; 
%i=find(TT>2*pi); TT(i)=TT(i)-2*pi;
%i=find(TT>2*pi); TT(i)=TT(i)-2*pi;

tt6=tt6-tt6(1);
i=find(tt6<0);   tt6(i)=tt6(i)+2*pi;
tt6(length(tt6))=2*pi;

%TT6=TT6-TT6(1);
%i=find(TT6<0);   TT6(i)=TT6(i)+2*pi;
%TT6(:,size(TT6,2))=2*pi;


if(pileheight>0)
    x6minmax=minmax(x6)+[-10 10];
    y6minmax=minmax(y6)+[-10 10];
    ix6=find((x6minmax(1)<x(:,1))&(x(:,1)<x6minmax(2)));
    iy6=find((y6minmax(1)<y(1,:))&(y(1,:)<y6minmax(2)));
    xxx=x(ix6,iy6); yyy=y(ix6,iy6); zzz=z(ix6,iy6);

    %z7=griddata(xxx,yyy,zzz,x6,y6);
    z7=linear_interp_from_grid(xxx,yyy,zzz,x6,y6);
    %z7=interpfromgrid(xxx,yyy,zzz,x6,y6);
    r8=rnd*rr6;
    t8=ones(size(rnd))*tt6;
    %x8=rnd*(x6-xpilecenter)+xpilecenter;
    %y8=rnd*(y6-ypilecenter)+ypilecenter;
    z8=rnd*(z7-zpilecenter)+zpilecenter;

    ir=(2:Nrnd-1);
    irm1=ir-1;
    irp1=ir+1;
    Nt=length(theta);
    it=1:Nt;
    itm1=[Nt-1 1:Nt-1];
    itp1=[2:Nt 2];
    Niter=max(Nrnd,Ntheta)*10;
    for iter=1:Niter %smoothing
        z8(ir,it)=0.25*(z8(irm1,it)+z8(irp1,it)+z8(ir,itm1)+z8(ir,itp1));
    end  
    h8=(1-rnd.^2)*ones(size(x6))*pileheight;
  
    [z9,h9]=linear_interp_from_grid(r8,t8,z8,rr,tt,h8);
    %z9=interpfromgrid(r8,t8,z8,rr,tt);
    %h9=interpfromgrid(r8,t8,h8,rr,tt);

    %z9=griddata(x8,y8,z8,x,y);
    %h9=griddata(x8,y8,h8,x,y);

    i=find(isfinite(z9));
    z(i)=z9(i); %+h9(i);
    h9(find(~isfinite(h9)))=0.0;
  
    %phdh=(1-rndct.^2-rndst.^2);
    %z6dz4=z6mmc(3)+rndct*dR6mm(3,1)+rndst*dR6mm(3,2);
    %surf(x,y,z); shading flat; axis image; colormap white; clh=camlight; hold on
    %plot3(x6,y6,z7,'-r');
    %view(0,90);
    %bob
end

Z7=griddata(X6,Y6,Z6,x,y);
Z7(find(isnan(Z7)))=-inf;
z=reshape(max(z(:),Z7(:)),size(x));

ix=find((380400<=x(:,1))&(x(:,1)<=381400));
iy=find((1846600<=y(1,:))&(y(1,:)<=1847600));

zz=z(ix,iy);

for iter=1:5000
  ifcontinue=0;
  for iiy=2:length(iy)-1
    iix=2:length(ix)-1;
    iixm1=iix-1;
    iixp1=iix+1;
    temp=[zz(iixm1,iiy) zz(iixp1,iiy) zz(iix,iiy-1) zz(iix,iiy+1)];
    temp2=minmax(temp);
    temp3=mean(temp,2);
    
    ii=find((zz(iix,iiy)<=temp2(:,1))&(zz(iix,iiy)<temp2(:,2)));
    if(~isempty(ii))
        ifcontinue=1;
    end
    zz(ii+1,iiy)=temp3(ii);
  end
  if(ifcontinue==0)
    break
  end
end

z(ix,iy)=zz;

if(ifplot)
    surf(x,y,z);
    %surf(x,y,z+p);
    shading flat; axis image; colormap white; clh=camlight; hold on
    if(pileheight>0)
        plot3(x6,y6,z7,'-r',x6,y6,z6,'-g');
    end
    hold off
    save debugme
    %axis([380000 382000 1846000 1848000]);
    %view(0,90);
    %view(-135,30)
    %bob
end
%clear zz Nz r theta Ntheta R Z1 THETA X1 Y2 rotang1 Z2 X2 rotang2 rotang3 X3 Y3 Z4 X4 Y4 Z5 i;

z=fixslope(x,y,z,65,10000);

return

function [zi,varargout]=linear_interp_from_grid(x,y,z,xi,yi,varargin)
    nout=max(nargout,1)-1;
    nin=max(nargin,5)-5;
    if((nin>0)&&(nout<nin))
        nin=nout;
    elseif(nout>nin+1)
        varargout(nin+2:nout)={[]};
    end
            
    %get triangles from rectangular grid "naturally"/quickly rather than 
    %from Delaunay Triangulation.
    [Nx,Ny]=size(x); %I have assumed that x, y, and z are the same size
    %more generally I will need to check this
    ix=shiftdim(repmat((1:Nx-1)',1   ,Ny-1),-1); %ix of grid "squares"
    iy=shiftdim(repmat((1:Ny-1) ,Nx-1,   1),-1); %iy of grid "squares"
    tri=zeros(6,Nx-1,Ny-1);
    tri(1,:,:)=(iy-1)*Nx+ix;
    tri(2,:,:)=tri(1,:,:)+1; %(iy-1)*Nx+ix+1;
    tri(3,:,:)=tri(2,:,:)+Nx;% iy   *Nx+ix+1;
    tri(4,:,:)=tri(3,:,:);
    tri(5,:,:)=tri(1,:,:);
    tri(6,:,:)=tri(1,:,:)+Nx; % iy   *Nx+ix;
    tri=reshape(tri,3,2*(Nx-1)*(Ny-1))';
    clear ix iy Nx Ny;
    %yay! I have the triangles

    %griddata would expect these to be vectors
    x=x(:);
    y=y(:);
    z=z(:);
    
    %griddata would expect these to be vectors
    siz=size(xi);
    xi=xi(:);
    yi=yi(:);
    
    %begin code copied from griddata
    
    % Find the nearest triangle (t)
    t = tsearch(x,y,tri,xi,yi);

    % Only keep the relevant triangles.
    out = find(isnan(t));
    if ~isempty(out), t(out) = ones(size(out)); end
    tri = tri(t,:);

    % Compute Barycentric coordinates (w).  P. 78 in Watson.
    del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
          (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
    w(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi) - ...
              (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
    w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi) - ...
              (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
    w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi) - ...
              (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
    w(out,:) = zeros(length(out),3);

    z = z(:).'; % Treat z as a row so that code below involving
                % z(tri) works even when tri is 1-by-3.
    zi = sum(z(tri) .* w,2);

    zi = reshape(zi,siz);

    if ~isempty(out), zi(out) = NaN; end

    if(nout==nin+1)
        varargout(nout)={out};
    end
    
    for iin=1:nin
       zz=varargin{iin};
       zz=zz(:).';
       varargout(iin)={reshape(sum(zz(tri) .* w,2),siz)};       
    end
    
    %save debugme
    %bob;
return

function z=fixslope(x,y,z,maxang,Maxloop)
%disp('entering fixslope.m');
NX=size(x,1);
NY=size(x,2);

IX=2:NX-1; IY=2:NY-1;
slope=zeros(NX,NY);
slope(IX,IY)=atan((...
    max(((z(IX+1,IY)-z(IX  ,IY))./(x(IX+1,IY)-x(IX  ,IY))).^2,...
        ((z(IX  ,IY)-z(IX-1,IY))./(x(IX  ,IY)-x(IX-1,IY))).^2)+...
    max(((z(IX,IY+1)-z(IX,IY  ))./(y(IX,IY+1)-y(IX,IY  ))).^2,...
        ((z(IX,IY  )-z(IX,IY-1))./(y(IX,IY  )-y(IX,IY-1))).^2)...
).^0.5)*180/pi;

%maxang=65;
[maxslope,imax]=max(slope(:));
if(maxslope>maxang)
  z2=z;
  iloop=0;
  while((maxslope>maxang)&&(iloop<Maxloop))
    iloop=iloop+1;
    %disp(sprintf('fixslope.m iloop=%g maxslope=%g at (%g,%g,%g) %g<=x<=%g %g<=y<=%g',iloop,maxslope,x(imax),y(imax),z(imax),minmax(x(:)'),minmax(y(:)')));
    for ix=2:NX-1
      iy=find(slope(ix,:)>maxang);
      if(length(iy))
        if(iy(1)==1) 
          iy=iy(2:length(iy));
        end
        if(iy(length(iy))==NY)
          iy=iy(1:length(iy)-1);
        end
      end
          
      z2(ix,iy)=(z(ix,iy+1)+...
                 z(ix-1,iy  )+z(ix+1,iy  )+...
                 z(ix,iy-1))/4;
    end
    z=z2;
    slope(IX,IY)=atan((...
        max(((z(IX+1,IY)-z(IX  ,IY))./(x(IX+1,IY)-x(IX  ,IY))).^2,...
            ((z(IX  ,IY)-z(IX-1,IY))./(x(IX  ,IY)-x(IX-1,IY))).^2)+...
        max(((z(IX,IY+1)-z(IX,IY  ))./(y(IX,IY+1)-y(IX,IY  ))).^2,...
            ((z(IX,IY  )-z(IX,IY-1))./(y(IX,IY  )-y(IX,IY-1))).^2)...
    ).^0.5)*180/pi;
    [maxslope,imax]=max(slope(:));
  end
end

%disp('exiting fixslope.m');

return;

function mm=minmax(a)
  mm=[min(a) max(a)];
return
