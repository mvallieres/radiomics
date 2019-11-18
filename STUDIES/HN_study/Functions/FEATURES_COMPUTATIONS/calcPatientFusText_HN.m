function [textures] = calcPatientFusText_HN(sDataPT,sDataCT,nameROI,CTweight_mat,scale_mat,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% function [textures] = calcPatientFusText_HN(sDataPT,sDataCT,nameROI,CTweight_mat,scale_mat,Ng_mat)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes FUSED texture features for all patients, for all 
% different combinations of the following texture extraction parameters:
% - CT weight: Weight given to MRI wavelet low-pass sub-bands in the 
%               PET/CT fusion process. 
% - Scale: Resolution at which the ROI is isotropically resampled.
% - Ng: Number of gray-levels in the quantization process. 
%
% Different extraction parameters are passed as arrays or cells in the
% function in order to test all possible combinations. This function is 
% used for FUSED scans specifically. See Ref. [1,2] and 'prepareVolume.m' 
% for more details.
% -------------------------------------------------------------------------
% REFERENCE:
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
% 1. sDataPT: File organizing the PT data to analyze.
% 2. sDataCT: File organizing the CT data to analyze.
% 3. nameROI: Cell of strings specifying the name of the ROI to use for texture analysis
%             --> {'GTV'}
% 4. CTweight_mat: Array vector specifying the different CT weights to test
%                  --> Ex: [1/4,1/3,1/2,2/3,3/4]
% 5. scale_mat: Array vector specifying the different 'Scale' values to test.
%               --> Ex: [1,2,3,4,5]
% 6. Ng_mat: Array vector specifying the different 'Ng' values to test.
%            --> Ex: [8,16,32,64]
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Structure with Global, GLCM, GLRLM, GLSZM and NGTDM textures
%             for all combinations of extraction parameters.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: March 2016
% -------------------------------------------------------------------------
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

textName = {{'Variance','Skewness','Kurtosis'}, ...
            {'Energy','Contrast','Entropy','Homogeneity','Correlation','SumAverage','Variance','Dissimilarity','AutoCorrelation'}, ...
            {'SRE','LRE','GLN','RLN','RP','LGRE','HGRE','SRLGE','SRHGE','LRLGE','LRHGE','GLV','RLV'}, ...
            {'SZE','LZE','GLN','ZSN','ZP','LGZE','HGZE','SZLGE','SZHGE','LZLGE','LZHGE','GLV','ZSV'}, ...
            {'Coarseness','Contrast','Busyness','Complexity','Strength'}};
        

% INITIALIZATION
textures.List = '';
textures.Parameters.CTweight = CTweight_mat;
textures.Parameters.Scale = scale_mat;
textures.Parameters.Algo = algo_cell;
textures.Parameters.Ng = Ng_mat;
nExperiment = numel(CTweight_mat)*numel(scale_mat)*numel(algo_cell)*numel(Ng_mat);
maskPT = zeros(size(sDataPT{2}.scan.volume.data)); maskCT = zeros(size(sDataCT{2}.scan.volume.data));
[contourVectPT] = findContours_HN(sDataPT,{nameROI},'PET'); [contourVectCT] = findContours_HN(sDataCT,{nameROI},'CT');
for c = 1:numel(contourVectPT)
    tempPT = zeros(size(sDataPT{2}.scan.volume.data)); tempCT = zeros(size(sDataCT{2}.scan.volume.data));
    [~,ROImaskPT,bbPT] = getROI(sDataPT,contourVectPT(c)); [~,ROImaskCT,bbCT] = getROI(sDataCT,contourVectCT(c));
    tempPT(bbPT(1,1):bbPT(1,2),bbPT(2,1):bbPT(2,2),bbPT(3,1):bbPT(3,2)) = ROImaskPT; tempCT(bbCT(1,1):bbCT(1,2),bbCT(2,1):bbCT(2,2),bbCT(3,1):bbCT(3,2)) = ROImaskCT;
    maskPT = maskPT + tempPT; maskCT = maskCT + tempCT;
end
maskPT(maskPT > 0) = 1; maskCT(maskCT > 0) = 1;
ROIbox_PT = sDataPT{2}.scan.volume.data; ROIbox_CT = sDataCT{2}.scan.volume.data;
bbPT = computeBoundingBox(maskPT); bbCT = computeBoundingBox(maskCT);
ROIbox_PT = ROIbox_PT(bbPT(1,1):bbPT(1,2),bbPT(2,1):bbPT(2,2),bbPT(3,1):bbPT(3,2)); ROIbox_CT = ROIbox_CT(bbCT(1,1):bbCT(1,2),bbCT(2,1):bbCT(2,2),bbCT(3,1):bbCT(3,2));
maskPT = maskPT(bbPT(1,1):bbPT(1,2),bbPT(2,1):bbPT(2,2),bbPT(3,1):bbPT(3,2));
szPT = size(ROIbox_PT); szCT = size(ROIbox_CT);
if szPT(3) ~= szCT(3)
    ROIbox_CT = imresize3D(ROIbox_CT,[],[szCT(1),szCT(2),szPT(3)],'cubic','fill');
end
pixelW = sDataPT{2}.scan.volume.spatialRef.PixelExtentInWorldX; sliceS = sDataPT{2}.scan.volume.spatialRef.PixelExtentInWorldZ;
scanType = 'Other';
wavelet = 'sym8'; % By default


experiment = 0;
fprintf('\n')
for w = 1:numel(CTweight_mat)
    [fusedBox] = fusePETCT(ROIbox_PT,ROIbox_CT,maskPT,'NoInv',CTweight_mat(w),wavelet);
    for s = 1:numel(scale_mat)
        try
            [ROIonly,~,ROIbox,maskBox] = prepareVolume(fusedBox,maskPT,scanType,pixelW,sliceS,1,scale_mat(s),'Global');
        end
        try
            [Global_text] = getGlobalTextures(ROIonly,100);
        catch
            for i = 1:numel(textName{1})
                Global_text.(textName{1}{i}) = 0;
            end
        end
        clear ROIonly
        for a = 1:numel(algo_cell)
            for n = 1:numel(Ng_mat)
                tic
                experiment = experiment + 1;
                strExperiment = ['Experiment',num2str(experiment)];
                nameExperiment = ['CTinv=','NoInv',', CTweight=',num2str(CTweight_mat(w),'%.2f'),', R=',num2str(1,'%.2f'),', Scale=',num2str(scale_mat(s)),', Quant.algo=',algo_cell{a},', Ng=',num2str(Ng_mat(n))];
                textures.List.(strExperiment) = nameExperiment;
                fprintf(['PERFORMING EXPERIMENT %u OF %u: ''',nameExperiment,''' ... '],experiment,nExperiment)
                try
                    [ROIonly,levels] = prepareVolume(ROIbox,maskBox,scanType,pixelW,pixelW,1,'pixelW','Matrix',algo_cell{a},Ng_mat(n)); % Pre-processing, WBPF and resampling already applied. Thus, we insert 'Other', R=1, sliceS = pixelW / Scale='pixelW', respectively
                end
                try
                    [GLCM] = getGLCM(ROIonly,levels); [GLCM_text] = getGLCMtextures(GLCM);
                catch
                    for i = 1:numel(textName{2})
                        GLCM_text.(textName{2}{i}) = 0;
                    end
                end
                try
                    [GLRLM] = getGLRLM(ROIonly,levels); [GLRLM_text] = getGLRLMtextures(GLRLM);
                catch
                    for i = 1:numel(textName{3})
                        GLRLM_text.(textName{3}{i}) = 0;
                    end            
                end
                try
                    [GLSZM] = getGLSZM(ROIonly,levels); [GLSZM_text] = getGLSZMtextures(GLSZM);
                catch
                    for i = 1:numel(textName{4})
                        GLSZM_text.(textName{4}{i}) = 0;
                    end                        
                end
                try
                    [NGTDM,countValid] = getNGTDM(ROIonly,levels); [NGTDM_text] = getNGTDMtextures(NGTDM,countValid);
                catch
                    for i = 1:numel(textName{5})
                        NGTDM_text.(textName{5}{i}) = 0;
                    end              
                end
                textures.(strExperiment).Global = Global_text;
                textures.(strExperiment).GLCM = GLCM_text;
                textures.(strExperiment).GLRLM = GLRLM_text;
                textures.(strExperiment).GLSZM = GLSZM_text;
                textures.(strExperiment).NGTDM = NGTDM_text;
                textures.(strExperiment).parameters.CTinv = 'NoInv';
                textures.(strExperiment).parameters.CTweight = CTweight_mat(w);
                textures.(strExperiment).parameters.R = 1;
                textures.(strExperiment).parameters.Scale = scale_mat(s);
                textures.(strExperiment).parameters.Algo = algo_cell{a};
                textures.(strExperiment).parameters.Ng = Ng_mat(n);
                clear ROIonly levels
                toc
            end
        end
        clear ROIbox maskBox
    end
    clear fusedBox
end

end