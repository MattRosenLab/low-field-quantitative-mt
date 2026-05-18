function [header] = Read_Tecmag_Header(filename)
header=[];
if (nargin<1)
    [filename,path]=uigetfile('*.tnt');
    filename=strcat(path,filename);
end
    
    if(exist(filename)==2)
    s=dir(filename);
    if s.bytes>1056
        Hdr_var=Read_Tecmag_hdr;

        Hdr_name=Hdr_var(:,1:28);
        Hdr_offset=Hdr_var(:,29:34);
        Hdr_type=Hdr_var(:,36:43);
        Hdr_size=Hdr_var(:,45:47);
        Hdr_desc=Hdr_var(:,48:85);
        % function strrep pour enlever les espaces.
    
        Number_var=size(Hdr_name);
        
        File_pt=fopen(filename,'r','l');
        for i=1:Number_var(1,1)
            Var_name=strrep(Hdr_name(i,:),' ','');
            Var_offset=eval(Hdr_offset(i,:));
            Var_type=strrep(Hdr_type(i,:),' ','');
            Var_size=eval(Hdr_size(i,:));
            Var_desc=Hdr_desc(i,:);
              
            fseek(File_pt,Var_offset,'bof');
            Var_data=fread(File_pt,Var_size,Var_type)';
            if (strcmp(Var_type,'char'))&((Var_size-1)>0)
                Var_data=char(Var_data);
            end
            eval(strcat('header.',Var_name,'=Var_data;'));
        end
        fclose(File_pt);
    end
end