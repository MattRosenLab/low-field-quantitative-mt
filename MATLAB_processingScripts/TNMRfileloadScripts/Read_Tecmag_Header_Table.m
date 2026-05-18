function arrays = Read_Tecmag_Header_Table(filename,Var_names)
% Modified from Read_Tecmag_Header.m by Dave Korenchan
header=[];
    
if exist(filename,'file')
    s=dir(filename);
    if s.bytes>1056
%         Hdr_var=Read_Tecmag_hdr;
% 
%         Hdr_name=Hdr_var(:,1:28);
%         Hdr_offset=Hdr_var(:,29:34);
%         Hdr_type=Hdr_var(:,36:43);
%         Hdr_size=Hdr_var(:,45:47);
%         Hdr_desc=Hdr_var(:,48:85);
%         % function strrep pour enlever les espaces.
%     
        Number_var=numel(Var_names);
        arrays=cell(Number_var,1);
        
        File_pt=fopen(filename,'r','l');
        remCarryOver=false;
        rem=[];
        while 1
            if remCarryOver
                line=rem;
                remCarryOver=false;
            else
                line=fgetl(File_pt);
            end
            for i=1:Number_var
                Var_name=Var_names{i};
                if strfind(line,Var_name)>0
                    %Make sure line is the beginning of the list! There
                    %should be a non-blank character after Var_name
                    rem=line(strfind(line,Var_name)+length(Var_name):end);
%                     if isempty(strfind(rem,[Var_name char(1)])) && ...
%                             isempty(strfind(rem,[char(0) Var_name char(1)]))
                    if ~strcmp(rem(1),{char(0),char(1)})
                        %Remove very end of rem -- that's the 1st entry
                        [~,val]=strtok(rem,char(0));
                        %ID whether there is a 'u', 'm', or 's' following
                        %the number. Keep timings in seconds by dividing
                        %appropriately
                        arrays{i}(1)=Parse_Tecmag_Variable(val);
                        ctr=1;
                        while 1
                            [lineval,rem]=strtok(fgetl(File_pt),char(0));
                            %Again, check for 'u', 'm', or 's' following
                            %the number. Keep in seconds
                            newval=Parse_Tecmag_Variable(lineval);                    
                            if ~isempty(newval) && ~isnan(newval)
                                ctr=ctr+1;
                                arrays{i}(ctr)=newval;
                            else
                                break
                            end
                            if ~isempty(rem) %use rem as next line, to see 
                                %if there's a variable of interest starting there!
                                remCarryOver=true;
                                break
                            end
                        end
%                         %Then remove variable from Var_names
%                         Var_names{i}=[];
                    end
                end
            end
            if ~ischar(line) || isempty(Var_names)
                break
            end
                
                
%                 Var_offset=eval(Hdr_offset(i,:));
%                 Var_type=strrep(Hdr_type(i,:),' ','');
%                 Var_size=eval(Hdr_size(i,:));
%                 Var_desc=Hdr_desc(i,:);
% 
%                 fseek(File_pt,Var_offset,'bof');
%                 Var_data=fread(File_pt,Var_size,Var_type)';
    %             if (strcmp(Var_type,'char'))&((Var_size-1)>0)
    %                 Var_data=char(Var_data);
    %             end
%                 eval(strcat('header.',Var_name,'=Var_data;'));
        end
        fclose(File_pt);
    end
end