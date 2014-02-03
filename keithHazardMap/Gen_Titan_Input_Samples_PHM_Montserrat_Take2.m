function Gen_Titan_Input_Samples_PHM_Montserrat_Take2(Nymacro)
    ifplot=1;

    if(ischar(Nymacro))
      Nymacro=str2num(Nymacro);
    end
    NymacroWant=Nymacro;
    Nymacro=2^round(log2(Nymacro));
    if(Nymacro~=NymacroWant)
        warning('Nymacro rounded from %g to %g',NymacroWant,Nymacro);
    end
    Ndiminmacro=4;
    
    NCellsPerEdge=2^ceil(log2(Nymacro)/Ndiminmacro);
    SafetyFactor=(NCellsPerEdge/2+1)/(NCellsPerEdge/2); %the purpose of the
    %SafetyFactor is to make the region we simulate in a little bigger than
    %the region we emulate in so that hopefully all/(most) of the resample
    %inputs will lie inside the convex hull of the points we
    %sampled/simulated.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MinEvent=5;
    MaxEvent=10; %9;
    %MiddleEvent=0.5*(MinEvent+MaxEvent);
    minlog10vol=MinEvent; %+SafetyFactor*(MinEvent-MiddleEvent);
    maxlog10vol=MaxEvent; %+SafetyFactor*(MaxEvent-MiddleEvent);
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
    %but for simulations cover the 0 to 360 range uniformly so the
    %simulations can be used for both angle draw distributions
    Mindir=0;
    Maxdir=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bedfriction angle uniformly distributed between 10 and 36 degrees
    BEDMIN=4.8;
    BEDMAX=12.2;
    %INT=BED+17+7*rand
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rand('twister',5489)
    %rand('seed',0);

%    fid=fopen('macro_resamples.tmp','w');
%    fprintf(fid,'additional file format lines=%g\n',3);
%
%    fprintf(fid,'%%Ndiminmacro=4: #of macro-emulator inputs (log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
%    fprintf(fid,'%%Nxmacro: the number of macro emulator input points to evaluate: [1] integer\n');
%    fprintf(fid,'%%{{x=(log10(volume [m^3]),Estart [UTME],Nstart [UTMN],BedFrictAng [deg]): [Nxmacro Ndiminmacro] doubles},{w: relative weight for hazmap assembly, sum(w)~=1 is ok: [1] double}}\n');


    r=BinOptLHSRand(Ndiminmacro,Nymacro);
    %r=rand(Ndiminmacro,Nxmacro)';
    
    log10vol=(maxlog10vol-minlog10vol)*r(:,1)+minlog10vol;

    direction=360*r(:,2);
    
     
    BEDFRICTANG=BEDMIN+(BEDMAX-BEDMIN)*r(:,3);
    INTFRICTANG=BEDFRICTANG+17+8*r(:,4);
    %maxbedfrict=max(BEDFRICTANG)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fid=fopen('uncertain_input_list.txt','w');
    fprintf(fid,'additional file format lines=3\n');
    fprintf(fid,'%%Ndiminmacro=4: of macro-emulator input dimensions (log10(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [1] integer\n');
    fprintf(fid,'%%Nymacro: number of simulations/mini-emulators: [1] integer\n');
    fprintf(fid,'%%y=(volume [m^3]),Direction [deg CC from east],BedFrictAng [deg],IntFrictAng [deg]): [Nymacro Ndiminmacro] array of doubles\n');
    fprintf(fid,'%g\n',Ndiminmacro,Nymacro);
    fprintf(fid,'%.10g %.10g %.10g %.10g\n',[log10vol direction BEDFRICTANG INTFRICTANG]');
    fclose(fid);
     
    if(0)
        fig=figure;
        pos=get(fig,'position');
        pos(2)=pos(2)-0.5*pos(4);
        pos(4)=pos(4)*1.5;
        set(fig,'position',pos,'paperunits','inches','paperposition',[.75 1 7 9]);
        
        subplot(3,2,1);
        plot(log10vol,STARTUTME,'bx'); axis square; axis([minlog10vol maxlog10vol STARTUTMECEN+STARTRADIUSMAX*[-1 1]]);
        xlabel('log10(vol)');
        ylabel('UTME');

        subplot(3,2,2);
        plot(log10vol,STARTUTMN,'bx'); axis square; axis([minlog10vol maxlog10vol STARTUTMNCEN+STARTRADIUSMAX*[-1 1]]);
        xlabel('log10(vol)');
        ylabel('UTMN');

        subplot(3,2,3);
        plot(log10vol,BEDFRICTANG,'bx'); axis square; axis([minlog10vol maxlog10vol BEDMIN BEDMAX]);
        xlabel('log10(vol)');
        ylabel('\phi_{bed} [deg]');
        
        subplot(3,2,4);
        plot(STARTUTME,STARTUTMN,'bx'); axis equal; axis([STARTUTMECEN+STARTRADIUSMAX*[-1 1] STARTUTMNCEN+STARTRADIUSMAX*[-1 1]]);
        xlabel('UTME');
        ylabel('UTMN');
        
        subplot(3,2,5);
        plot(STARTUTME,BEDFRICTANG,'bx'); axis square; axis([STARTUTMECEN+STARTRADIUSMAX*[-1 1] BEDMIN BEDMAX]);
        xlabel('UTME');
        ylabel('\phi_{bed} [deg]');
                
        subplot(3,2,6);
        plot(STARTUTMN,BEDFRICTANG,'bx'); axis square; axis([STARTUTMNCEN+STARTRADIUSMAX*[-1 1] BEDMIN BEDMAX]);
        xlabel('UTMN');
        ylabel('\phi_{bed} [deg]');        
    end
        
return;