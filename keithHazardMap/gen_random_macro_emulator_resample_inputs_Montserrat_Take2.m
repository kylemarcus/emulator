function gen_random_macro_emulator_resample_inputs_Montserrat_Take2(Nxmacro)
    ifplot=0;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %these 2 parameters are from the SAMSI technometrics paper for
    %Montserrat
    lambda=12425; %events per year
    alpha=0.65; %pareto exponent, unitless    
    Years=10;
    minlog10vol=5;
    maxlog10vol=log10(3*10^9);
    
    %probability of log10volume=
    %Years*lambda*alpha*log(10)*(vol.^-alpha).*exp(-Years*lambda*vol.^-alpha);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %direction historical frequency
    histangle=[...
         20 573;...
         78 151;...
         90  17;...
         98 230;...
        125  37;...
        143  87;...
        180  89;...
        230   1;...
        270  80];
    NHist=sum(histangle(:,2));
    SumHistAng=cumsum(histangle(:,2));
    IHistAng=repmat(length(histangle),NHist,1);
    for ii=length(histangle)-1:-1:1
        IHistAng(1:SumHistAng(ii))=ii;
    end
    rand('twister',5489);
    IHistAng=IHistAng(randperm(NHist));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bedfriction angle uniformly distributed between 5 and 12 degrees
    BEDMIN=5;
    BEDMAX=12;
    %INT=BED+17+7*rand
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ischar(Nxmacro))
        Nxmacro=str2num(Nxmacro);
    end
    NxmacroWant=Nxmacro;
    Nxmacro=2^round(log2(Nxmacro));
    if(Nxmacro~=NxmacroWant)
        warning('Nxmacro rounded from %g to %g',NxmacroWant,Nxmacro);
    end
    
    rand('twister',5489)
    %rand('seed',0);

    fid=fopen('macro_resamples.tmp','w');
    fprintf(fid,'additional file format lines=%g\n',3);

    fprintf(fid,'%%Ndiminmacro=4: #of macro-emulator inputs (log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
    fprintf(fid,'%%Nxmacro: the number of macro emulator input points to evaluate: [1] integer\n');
    fprintf(fid,'%%{{x=(log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [Nxmacro Ndiminmacro] doubles},{w: relative weight for hazmap assembly, sum(w)~=1 is ok: [1] double}}\n');

    Ndiminmacro=4;
    r=BinOptLHSRand(Ndiminmacro,Nxmacro);
    %r=rand(Ndiminmacro,Nxmacro)';
   
    log10vol=(maxlog10vol-minlog10vol)*r(:,1)+minlog10vol;
    vol=10.^log10vol;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %probability of log10(volume)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    pdraw=1/(maxlog10vol-minlog10vol); 
    pfrechet=Years*lambda*alpha*log(10)*(vol.^-alpha).*exp(-Years*lambda*vol.^-alpha);  


    w=pfrechet./pdraw; %non normalized likelihood ratio for importance sampling
    sumwdivNxmacro=sum(w)/Nxmacro
    w=w/mean(w);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save DEBUGME;
    Direction=histangle(IHistAng(ceil(NHist*r(:,2))),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Friction Angles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BEDFRICTANG=BEDMIN+(BEDMAX-BEDMIN)*r(:,3);
    INTFRICTANG=BEDFRICTANG+17+7*r(:,4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(fid,'%d\n%d\n',Ndiminmacro,Nxmacro);
    fprintf(fid,'%.10g %.10g %.10g %.10g %.10g\n',[log10vol Direction BEDFRICTANG INTFRICTANG w]');
    fclose(fid);
    
    if(ifplot)
        fig=figure;
        pos=get(fig,'position');
        pos(2)=pos(2)-0.5*pos(4);
        pos(4)=pos(4)*1.5;
        set(fig,'position',pos,'paperunits','inches','paperposition',[.75 1 7 9]);
        
        subplot(3,2,1);
        plot(log10vol,Direction,'bx'); axis square; %axis([minlog10vol maxlog10vol STARTUTMECEN+STARTRADIUSMAX*[-1 1]]);
        xlabel('log10(vol)');
        ylabel('Direction [deg]');

        subplot(3,2,2);
        plot(log10vol,BEDFRICTANG,'bx'); axis image; %square; axis([minlog10vol maxlog10vol BEDMIN BEDMAX]);
        xlabel('log10(vol)');
        ylabel('\phi_{bed} [deg]');

        subplot(3,2,3);
        plot(log10vol,INTFRICTANG,'bx'); axis image; %square; axis([minlog10vol maxlog10vol BEDMIN BEDMAX]);
        xlabel('log10(vol)');
        ylabel('\phi_{int} [deg]'); 

        subplot(3,2,4);
        plot(Direction,BEDFRICTANG,'bx'); axis square; %axis([STARTUTMECEN+STARTRADIUSMAX*[-1 1] BEDMIN BEDMAX]);
        xlabel('Direction [deg]');
        ylabel('\phi_{bed} [deg]');        
        
        subplot(3,2,5);
        plot(Direction,INTFRICTANG,'bx'); axis square; %axis([STARTUTMECEN+STARTRADIUSMAX*[-1 1] BEDMIN BEDMAX]);
        xlabel('Direction');
        ylabel('\phi_{int} [deg]');
                
        subplot(3,2,6);
        plot(BEDFRICTANG,INTFRICTANG,'bx'); axis image; %square; axis([STARTUTMNCEN+STARTRADIUSMAX*[-1 1] BEDMIN BEDMAX]);
        xlabel('\phi_{bed} [deg]');
        ylabel('\phi_{int} [deg]');        
    end
        
return;