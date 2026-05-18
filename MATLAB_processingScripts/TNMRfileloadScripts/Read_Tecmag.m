function [Ms,Header,Var_data,sw]=Read_Tecmag(filename,dsize)
    Ms=[];
    Header=[];

    if (nargin<1)
%         [filename,path]=uigetfile();
        [filename,path]=uigetfile('*.tnt');
        filename=strcat(path,filename);
        
    end
     Header=Read_Tecmag_Header(filename);
     ds=Header.acq_points;
     ss=Header.actual_npts;
    if (nargin<2)
        nc=ss(1)/ds(1);
        dsize=[nc ds(1) ss(2:end)];
    end
    
     Hdr_var=Read_Tecmag_hdr;

        Hdr_name=Hdr_var(:,1:28);
        Hdr_offset=Hdr_var(:,29:34);
        Hdr_type=Hdr_var(:,36:43);
        Hdr_size=Hdr_var(:,45:47);
        Hdr_desc=Hdr_var(:,48:85);
        % function strrep pour enlever les espaces.
    
     Number_var=size(Hdr_name);
        
        File_pt=fopen(filename,'r','l');
        for i=20
            Var_name=strrep(Hdr_name(i,:),' ','');
            Var_offset=eval(Hdr_offset(i,:));
            Var_type=strrep(Hdr_type(i,:),' ','');
            Var_size=eval(Hdr_size(i,:));
            Var_desc=Hdr_desc(i,:);
              
            fseek(File_pt,Var_offset,'bof');
            Var_data=fread(File_pt,Var_size,Var_type)';
        end  
    
    
% lecture des données
offset=1056;
   Ms=Read_Raw_Field(filename,offset,[2,dsize],'float','l');
   Ms=Ms(1,:,:,:,:)+1i*Ms(2,:,:,:,:);
   Ms=squeeze(Ms);
   %Ms=permute(Ms,[2 1 3 4]);
   
   fclose(File_pt);
   
   

   
   
