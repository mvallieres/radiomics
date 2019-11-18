function organizeSeparateTextures_HN(pathWORK,nPatient,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% function organizeSeparateTextures_HN(pathWORK,nPatient,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes SEPARATE texture features for all patients, for
% all different combinations of the following texture extraction parameters:
% - WPBF ratio 'R': Ratio specifying the ratio of weight given to band-pass 
%                   coefficients to the weigth given to the rest of 
%                   coefficients (HHH and LLL) in the wavelet band-pass 
%                   filtering process.
% - Scale: Resolution at which the ROI is isotropically resampled.
% - Quant. Algo.: Quantization algorithm used on analyzed ROI.
% - Ng: Number of gray-levels in the quantization process. 
%
% Different extraction parameters are passed as arrays or cells in the
% function in order to read all possible combinations. This function is 
% used for SEPARATE scans specifically. See Ref. [1,2] and 'prepareVolume.m' 
% for more details.
%
% Textures are organized in two cells: 'textures_PET'and textures_CT, where
% each cell entry organizes the data for all patients for a specific 
% extraction parameter combination.
%
% -- > For example, size(textures_PET) is equal to:
% [numel(R_mat),numel(scale_cell),numel(algo_cell),numel(Ng_mat)]
%
% In addition, the Spearman's rank correlation between each texture 
% features and each different outcome is computed and saved in the cells 
% for each texture feature of each entry. Results are saved in the HN 
% WORKSPACE with the names given above.
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Vallieres, M. et al. (2015). FDG-PET/CT radiomics models for the 
%     early prediction of different tumour outcomes in head and neck cancer.
%     The Journal of Nuclear Medicine, aa(bb), xxx-yyy. 
%     doi:
% [2] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - pathWORK: Full path to the HN WORKSPACE directory.
% - nPatient: Number of patients to read.
% - R_mat: Array vector specifying the different 'R' ratio values to test.
% - scale_cell: Cell vector specifying the different 'Scale' values to test.
% - algo_cell: Cell vector specifying the different 'Quant. Algo.' values to test.
% - Ng_mat: Array vector specifying the different 'Ng' values to test.
% -------------------------------------------------------------------------
% EXAMPLE:
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Uniform','Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2015
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
cd(pathWORK), load('outcomes')

% INITIALIZATION
scans = {'PET','CT'};
nScans = numel(scans);
suffix = ['_TEXT'];
nameOutcomes = fieldnames(outcomes);
nOutcomes = length(nameOutcomes);

for scan = 1:numel(scans)
    cd([pathWORK,'/TEXTURES'])
    temp = load('Patient1_PET_TEXT'); temp = struct2cell(temp); temp = temp{1}; % In order to get the necessary 'nameType' and 'nameFeature' fields
    nameType = fieldnames(temp.Experiment1); nType = numel(nameType); % All texture types are the same
    text = cell(numel(R_mat),numel(scale_cell),numel(algo_cell),numel(Ng_mat)); % Initialization of master cell
    tempText = cell(1,nPatient); % Cell used to load patient textures only once
    for p = 1:nPatient
        load(['Patient',num2str(p),'_',scans{scan},suffix]) % Variable 'textures' is now in MATLAB workspace
        tempText{p} = textures;
    end
    experiment = 0;
    for r = 1:numel(R_mat)
        for s = 1:numel(scale_cell)
            for a = 1:numel(algo_cell)
                for n = 1:numel(Ng_mat)
                    text{r,s,a,n} = struct;
                    experiment = experiment + 1;
                    strExperiment = ['Experiment',num2str(experiment)];
                    for t = 1:nType-1 % Last field is 'parameters'
                        nameFeature = fieldnames(temp.(strExperiment).(nameType{t})); nFeature = numel(nameFeature);
                        for f = 1:nFeature
                            data = zeros(nPatient,1);
                            for p = 1:nPatient
                                data(p,1) = tempText{p}.(strExperiment).(nameType{t}).(nameFeature{f});
                            end
                            text{r,s,a,n}.(nameType{t}).(nameFeature{f}).Data = data;
                            for o = 1:nOutcomes
                                text{r,s,a,n}.(nameType{t}).(nameFeature{f}).Spearman.(nameOutcomes{o}) = corr(data,outcomes.(nameOutcomes{o}),'type','Spearman');
                            end
                        end
                    end
                end
            end
        end
    end
    cd(pathWORK), save(['textures_',scans{scan}],'text')
end

cd(startpath)
end