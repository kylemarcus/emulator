function gen_random_macro_emulator_resample_inputs_Montserrat_Take2(Nxmacro)
    
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

    fid=fopen('macro_resamples.tmp','w');

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
    Direction=histangle(IHistAng(ceil(NHist*r(:,2))),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Friction Angles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BEDFRICTANG=BEDMIN+(BEDMAX-BEDMIN)*r(:,3);
    INTFRICTANG=BEDFRICTANG+17+7*r(:,4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N=length(log10vol,1);
    ind=1:1:N;
    fprintf(fid,'%d %.10g %.10g %.10g %.10g %.10g\n',[ind log10vol Direction BEDFRICTANG INTFRICTANG w]');

    fclose(fid);
    
return;
