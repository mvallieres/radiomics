function organizeFusedTextures_STS(pathWORK,nPatient,MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% function organizeFusedTextures(pathWORK,nPatient,MRIinv_cell,MRIweight_mat,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function organizes FUSED texture features computed for all 
% patients, for all different combinations of the following texture 
% extraction parameters:
% - MRI Inv.: Inversion of MRI intensities in the PET/MRI fusion process.
% - MRI weight: Weight given to MRI wavelet sub-bands in the PET/MRI fusion 
%               process.
% - WPBF ratio 'R': Ratio specifying the ratio of weight given to band-pass 
%                   coefficients to the weigth given to the rest of 
%                   coefficients (HHH and LLL) in the wavelet band-pass 
%                   filtering process.
% - Scale: Resolution at which 'volume' is isotropically resampled.
% - Quant. Algo.: Quantization algorithm used on 'volume'.
% - Ng: Number of gray-levels on the quantization process. 
%
% Different extraction parameters are passed as arrays or cells in the
% function in order to read all possible combinations. This function is 
% used for FUSED scans specifically. See Ref. [1] and 'prepareVolume.m' 
% for more details.
%
% Textures are organized in three cells: 'textures_PET_T1' and 
% 'textures_PET_T2FS', where each cell entry organizes the data for all 
% patients for a specific extraction parameter combination. T2-weighted 
% fat-saturated and STIR scans are combined in the same texture category 
% and simply refered to as 'T2FS' scans (T2-weighted fat-suppression scans).
% 
% -- > For example, size(textures_PET_T1) is equal to:
% [numel(MRIinv_cell),numel(MRIweight_mat),numel(R_mat),numel(scale_cell),numel(algo_cell),numel(Ng_mat)]
%
% In addition, the Spearman's rank correlation between each texture 
% features and lung metastases is computed and saved in the cells for each 
% texture feature of each entry. Results are saved in the STS WORKSPACE
% with the names given above.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the STS WORKSPACE directory.
% - nPatient: Number of patients to read.
% - MRIinv_cell: Cell vector specifying the different 'MRI Inv.' values to test.
% - MRIweight_mat: Array vector specifying the different 'MRI weight' values to test.
% - R_mat: Array vector specifying the different 'R' ratio values to test.
% - scale_cell: Cell vector specifying the different 'Scale' values to test.
% - algo_cell: Cell vector specifying the different 'Quant. Algo.' values to test.
% - Ng_mat: Array vector specifying the different 'Ng' values to test.
% -------------------------------------------------------------------------
% EXAMPLE:
% MRIinv_cell = {'NoInv','Inv'};
% MRIweight_mat = [1/4,1/3,1/2,2/3,3/4];
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres
%
%    This package is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this package.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

startpath = pwd;
cd(pathWORK), load('outcome')

% INITIALIZATION
scans = {'PET_T1','PET_T2FS'};
nScans = numel(scans);
suffix = ['_TEXT'];

fprintf('ORGANIZING TEXTURES FROM FUSED SCANS ... ')
for scan = 1:numel(scans)
    cd([pathWORK,'/TEXTURES'])
    temp = load('Patient1_PET_T1_TEXT'); temp = struct2cell(temp); temp = temp{1}; % In order to get the necessary 'nameType' and 'nameFeature' fields
    nameType = fieldnames(temp.Experiment1); nType = numel(nameType); % All texture types are the same
    text = cell(numel(MRIinv_cell),numel(MRIweight_mat),numel(R_mat),numel(scale_cell),numel(algo_cell),numel(Ng_mat)); % Initialization of master cell
    tempText = cell(1,nPatient); % Cell used to load patient textures only once
    for p = 1:nPatient
        load(['Patient',num2str(p),'_',scans{scan},suffix]) % Variable 'textures' is now in MATLAB workspace
        tempText{p} = textures;
    end
    experiment = 0;
    for i = 1:numel(MRIinv_cell)
        for w = 1:numel(MRIweight_mat)
            for r = 1:numel(R_mat)
                for s = 1:numel(scale_cell)
                    for a = 1:numel(algo_cell)
                        for n = 1:numel(Ng_mat)
                            text{i,w,r,s,a,n} = struct;
                            experiment = experiment + 1;
                            strExperiment = ['Experiment',num2str(experiment)];
                            for t = 1:nType-1 % Last field is 'parameters'
                                nameFeature = fieldnames(temp.(strExperiment).(nameType{t})); nFeature = numel(nameFeature);
                                for f = 1:nFeature
                                    data = zeros(nPatient,1);
                                    for p = 1:nPatient
                                        data(p,1) = tempText{p}.(strExperiment).(nameType{t}).(nameFeature{f});
                                    end
                                    text{i,w,r,s,a,n}.(nameType{t}).(nameFeature{f}).Data = data;
                                    text{i,w,r,s,a,n}.(nameType{t}).(nameFeature{f}).Spearman = corr(data,outcome,'type','Spearman');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    cd(pathWORK), save(['textures_',scans{scan}],'text')
end
fprintf('DONE\n')

cd(startpath)
end