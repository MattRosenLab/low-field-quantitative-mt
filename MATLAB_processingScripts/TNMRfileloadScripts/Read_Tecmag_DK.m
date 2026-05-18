function [Ms,Header,Var_data,filename]=Read_Tecmag_DK(filename,dsize,varargin)
% DK edit, 11/19/24: Return filename
Ms=[];
    Header=[];

    if (nargin<1)
        [filename,path]=uigetfile('/Users/mat/Documents/Matlab_MS/Mathieu/TNMR/DATA/*.tnt');
        filename=strcat(path,filename);
        
    end
    Header=Read_Tecmag_Header(filename);
    ds=Header.acq_points;
    ss=Header.actual_npts;
    if nargin < 3
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

            fclose(File_pt);  %added by B.Armstrong
        % lecture des données
        offset=1056;
           Ms=Read_Raw_Field(filename,offset,[2,dsize],'float','l');
           Ms=Ms(1,:,:,:,:)+1i*Ms(2,:,:,:,:);
           Ms=squeeze(Ms);
           %Ms=permute(Ms,[2 1 3 4]);

    else
        npnts=ss(2);
        % Determine channels
        channels = Read_Tecmag_channels(filename,Header.data_bytes);
        nc = length(channels);
        if nc == 1
            fprintf('Acquisition detected on channel %g\n',channels);
        else
            fprintf('Acquisition detected on channels %s\n',num2str(reshape(channels,1,nc)));
        end
        if (nargin<2)
    %         if max(Header.npts./Header.actual_npts) == 1
                %data has not been separated by Tnmr
                sorted_flag = false;

    %             nc=ss(1)/ds(1);
                dsize=[nc ds(1) ss(2:end)];
    %         else
    %             % data has been separated by Tnmr already
    %             sorted_flag = true;
    %             
    % %             nc = max(Header.npts./Header.actual_npts);
    %             dsize=[1 ds(1) ss(2:end)];
    %         end
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
    %             eval(['Var_data.' Var_name '=fread(File_pt,Var_size,Varr_type)''']);
            end  
        % Added by C. LaPierre March 29, 2012
        fclose(File_pt);    

        % lecture des données
        offset=1056;
           Ms=Read_Raw_Field(filename,offset,dsize,'float','l',varargin);

        % Sort data according to specified dimensions
        n_SA1 = length(find(channels <=4));
        n_SA2 = length(find(channels >=5));
        if sorted_flag
            Ms = reshape(Ms,ds(1),nc,npnts);
            Ms = permute(Ms,[1,3,2]);
        else
    %         switch nc
    %             case {1,2,3,4}
    %                 Ms2 = reshape(Ms,nc,ds(1),npnts);
    %                 Ms = permute(Ms2,[2,3,1]);
    %             case {5,6,7}
    %                 Ms2=zeros(nc,ds(1),npnts);
    %                 for loop = 1:npnts
    %                     Ms2(1:4,:,loop) = reshape(Ms(((loop-1)*nc*ds(1)+1):((loop-1)*nc*ds(1)+4*ds(1))),4,ds(1),1);
    %                     Ms2(5:nc,:,loop) = reshape(Ms(((loop-1)*nc*ds(1)+4*ds(1)+1):loop*nc*ds(1)),nc-4,ds(1),1);
    %                 end
    %                 Ms = permute(Ms2,[2,3,1]);
    %             case {8}
    %                 Ms2 = reshape(Ms,4,2*ds(1),npnts);
    % 
    %                 Ms3 = permute(Ms2,[2,3,1]);
    %                 Ms=cat(3,Ms3(1:end/2,:,:),Ms3(end/2+1:end,:,:));
    %         end
            Ms2=zeros(nc,ds(1),npnts);
            for loop = 1:npnts
                Ms2(1:n_SA1,:,loop) = reshape(Ms(((loop-1)*nc*ds(1)+1):((loop-1)*nc*ds(1)+n_SA1*ds(1))),n_SA1,ds(1),1);
                Ms2(n_SA1+1:nc,:,loop) = reshape(Ms(((loop-1)*nc*ds(1)+n_SA1*ds(1)+1):loop*nc*ds(1)),n_SA2,ds(1),1);
            end
            Ms = permute(Ms2,[2,3,1]);
        end
    end % end if 
end