 function [Z,offsets_kHz,Msat]=ZspecProc(Msat,B1_Hz) 

% 12/19/2024

% This function takes z-spectroscopy data and the corresponding saturation 
% B1 amplitudes, calculates the z-spectra, and displays them.
%
% INPUTS: 
%   Msat        --  Vector of water peak integrals obtained from running 
%                   NMR_UI_beta. Length is the total number of FIDs 
%                   ( = (# of sat offsets) x (# of sat powers) )   
%   B1_Hz       --  Vector of B1 amplitudes pertaining to the 3D table of
%                   F1 attenuation values, in units of Hz 
%
% OUTPUTS:
%   Z           --  Array of calculated z-spectra. 1st dimension is the
%                   total # of saturation offsets, ordered from lowest to 
%                   highest. 2nd dimension is the # of saturation 
%                   amplitudes. 
%   offsets_kHz --  Vector containing the loaded saturation offsets,
%                   ordered from lowest to highest.
%   Msat        --  Array of reordered integrals from input Msat. Same
%                   dimensions as Z.
%

% Have user specify the data file to determine the offsets used
[~,~,~,fname]=Read_Tecmag_DK;
disp('Reading in offset and power information from tables in header...')
tableVals=Read_Tecmag_Header_Table(fname,{'o1_0:2'});%,'at1:3'});
offsetsHz=tableVals{1}';
% dBattn=tableVals{2}';

% Reshape input data based upon the sizes of the offsets and powers read
% in, then reorder so that offsets are monotonically increasing
Msat=reshape(Msat,numel(offsetsHz),numel(B1_Hz));
[offsetsHz,sortIdx]=sort(offsetsHz,'ascend');
Msat=Msat(sortIdx,:);

% Calculate the z-spectrum for each power using the highest offset, in
% absolute value
% M0idx=find(abs(offsetsHz)==max(abs(offsetsHz)));
Z=Msat./repmat(Msat(abs(offsetsHz)==max(abs(offsetsHz)),:),size(Msat,1),1);

% Z(M0idx,:)=[];
% offsetsHz(M0idx)=[];

% Find positive and negative offset pairs for calculating MTR asymmetry
offsPosIdx=find(offsetsHz>0);
offsNegIdx=find(offsetsHz<0);
for ii=1:numel(offsPosIdx)
    offPos=offsetsHz(offsPosIdx(ii));
    if sum(find(abs(offsetsHz(offsPosIdx(ii))+offsetsHz(offsNegIdx))<1e-3,1))~=0
        offsPairsPosIdx(ii)=offsPosIdx(ii);
        offsPairsNegIdx(ii)=find(abs(offPos+offsetsHz(offsNegIdx))<1e-3);
    end
end

% Calculate MTR asymmetry
MTRasym=Z(offsPairsNegIdx,:)-Z(offsPairsPosIdx,:);

% Plot results
offsets_kHz=offsetsHz/1000;
offsMTR=offsets_kHz(offsPairsPosIdx);
leglbls=cell(size(Z,2),1);
for ii=1:size(Z,2)
    leglbls{ii}=['B_1 = ' num2str(B1_Hz(ii)) ' Hz'];
end
plotIdx=(abs(offsetsHz)~=max(abs(offsetsHz))); %plot all but M0

figure; plot(offsets_kHz(plotIdx),Z(plotIdx,:),'-*')
ylim([0 1.05])
xlabel('Offset (kHz)')
ylabel('M_{sat}/M_0')
legend(leglbls,'Location','southeast')
set(gca,'Xdir','reverse')
title('Z-spectra')

figure; plot(offsMTR,MTRasym,'-*')
xlabel('Offset (kHz)')
ylabel('MTR_{asym}')
legend(leglbls,'Location','southeast')
set(gca,'Xdir','reverse')
title('MTR asymmetry')
end