function ExtractSampleOverRideFile(samplenumber)
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
    y=fscanf(fid,'%g',[Ndimin Ny])';
    fclose(fid);
    make_montserrat_take2_pwem_dem(10^y(samplenumber,1),y(samplenumber,2),y(samplenumber,3),y(samplenumber,4),samplenumber);
    
    
    %fid=fopen(sprintf('VolENstartBedSample.%06g',samplenumber),'w');
    %fprintf(fid,'Nsample=%g\n',Ny);
    %fprintf(fid,'isample=%g\n',samplenumber);
    %fprintf(fid,'log10(Vol [m^3])=%.10g\nUTM East  Center=%.10g\nUTM North Center=%.10g\nBed Frict [deg] =%.10g',y(samplenumber,:));
    %fclose(fid);
return;
