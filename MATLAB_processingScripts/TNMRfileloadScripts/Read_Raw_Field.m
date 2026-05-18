function [M]=Read_Raw_Field(filename,Offset,Mat_size,Mat_type,Big_Little)
ok1=1;
ok2=1;

if nargin==0
filename=uigetfile;
end

if nargin<2
    % matrix dimension
    s=inputdlg({'Enter offset :','Enter Matrix dimension :'},...
                     'Input for Read_Raw_Field',...
                     1,...
                     {'0','[1]'});
    Offset=eval(char(s(1)));
    Mat_size=eval(char(s(2)));
    
    % type
    type_list={'schar','uchar','int8','int16','int32','int64','uint8','uint16','uint32','uint64','float32','float64','double',...
                         'char','short','int','long','ushort','uint','ulong','float',...
                         'bitN','ubitN'};
    [s,ok1]=listdlg('PromptString','Select a type:',...
                'SelectionMode','single',...
                'ListString',type_list);
    Mat_type=char(type_list(s));

    % number of bits if 'bitN'
        if (strcmp(Mat_type,'bitN')==1)|(strcmp(Mat_type,'ubitN')==1)
            s=inputdlg({'Enter number of bit(s):'},...
                'Read_Raw_Field Bit',...
                     1,...
                     {'1'});
            if  (strcmp(Mat_type,'bitN')==1) 
                Mat_type=['bit' char(s)];
            else
                Mat_type=['ubit' char(s)];
            end
        end
    
    % Endian mode
    Big_list={'cray','ieee-be','ieee-le','ieee-be.l64','ieee-le.l64','native','vaxd','vaxg'},
    [s,ok2]=listdlg('PromptString','Select a Big/Little Endian:',...
                'SelectionMode','single',...
                'ListString',Big_list);   
    Big_Little=char(Big_list(s));
end

M=[];

if ok1&ok2
    Msize=size(Mat_size);
    tot_size=1;
    for i=1:Msize(1,2)
        tot_size=tot_size*Mat_size(1,i);
    end

    fp=fopen(filename,'r',Big_Little);
    if fp~=-1
        fseek(fp,Offset,'bof');
        dat=fread(fp,tot_size,Mat_type);
        fclose(fp);
        M=dat;
        if (size(Mat_size)~=[1 1])|(Mat_size(1,1)>1)
            M=reshape(M,Mat_size);
        end  
    end
end