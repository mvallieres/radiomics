function [textures] = calcPatientFusText_HNissam(sDataPET,sDataCT,roiNumb,CTinv_cell,CTweight_mat,R_mat,scale_cell,algo_cell,Ng_mat,fusWeights)
% -------------------------------------------------------------------------
% function [textures] = calcPatientFusText_HN(sDataPET,sDataCT,roiNumb,CTinv_cell,CTweight_mat,R_mat,scale_cell,algo_cell,Ng_mat,weights)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes FUSED texture features for all patients, for all 
% different combinations of the following texture extraction parameters:
% - CT Inv.: Inversion of CT intensities in the PET/CT fusion process.
% - CT weight: Weight given to MRI wavelet low-pass sub-bands in the 
%               PET/CT fusion process. 
% - WPBF ratio 'R': Ratio specifying the ratio of weight given to band-pass 
%                   coefficients to the weigth given to the rest of 
%                   coefficients (HHH and LLL) in the wavelet band-pass 
%                   filtering process.
% - Scale: Resolution at which the ROI is isotropically resampled.
% - Quant. Algo.: Quantization algorithm used on analyzed ROI.
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
% - sDataPET: File organizing the PET data.
% - sDataCT: File organizing the CT data.
% - roiNumb: Array of size [N X 2] specifying the contour number 
%            corresponding to the GTV in sData{2}.scan.contour for all 
%            patients. First column is for the PET data, second column for 
%            CT data.
% - CTinv_cell: Cell vector specifying the different 'CT Inv.' values to test.
% - CTweight_mat: Array vector specifying the different 'CT weight' values to test.
% - R_mat: Array vector specifying the different 'R' ratio values to test.
% - scale_cell: Cell vector specifying the different 'Scale' values to test.
% - algo_cell: Cell vector specifying the different 'Quant. Algo.' values to test.
% - Ng_mat: Array vector specifying the different 'Ng' values to test.
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Structure with Global, GLCM, GLRLM, GLSZM and NGTDM textures
%             for all combinations of extraction parameters.
% -------------------------------------------------------------------------
% EXAMPLE:
% MRIinv_cell = {'NoInv','Inv'};
% MRIweight_mat = [1/4,1/3,1/2,2/3,3/4];
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
% [textures] = calcPatientFusText_STS(sDataPET,sDataCT,1,CTinv_cell,CTweight_mat,R_mat,scale_cell,algo_cell,Ng_mat);
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2015
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

% INITIALIZATION
textures.List = '';
textures.Parameters.CTinv = CTinv_cell;
textures.Parameters.CTweight = CTweight_mat;
textures.Parameters.R = R_mat;
textures.Parameters.Scale = scale_cell;
textures.Parameters.Algo = algo_cell;
textures.Parameters.Ng = Ng_mat;
nExperiment = numel(CTinv_cell)*numel(CTweight_mat)*numel(R_mat)*numel(scale_cell)*numel(algo_cell)*numel(Ng_mat);
ROIbox_PET = getROIbox(sDataPET,roiNumb); ROIbox_CT = getROIbox(sDataCT,roiNumb);
mask = sDataPET{2}.scan.contour(roiNumb).boxMask;
pixelW = sDataPET{2}.scan.pixelW; sliceS = sDataPET{2}.scan.sliceS; scanType = 'Other';
wavelet = 'sym8'; % By default


experiment = 0;
fprintf('\n')
for i = 1:numel(CTinv_cell)
    for w = 1:numel(CTweight_mat)
        [fusedBox] = fusePETCTissam(ROIbox_PET,ROIbox_CT,mask,CTinv_cell{i},CTweight_mat(w),wavelet,fusWeights);
        for r = 1:numel(R_mat)
            for s = 1:numel(scale_cell)
                [ROIonly,~,ROIbox,maskBox] = prepareVolume(fusedBox,mask,scanType,pixelW,sliceS,R_mat(r),scale_cell{s},'Global');
                [Global_text] = getGlobalTextures(ROIonly,100);
                for a = 1:numel(algo_cell)
                    for n = 1:numel(Ng_mat)
                        tic
                        experiment = experiment + 1;
                        strExperiment = ['Experiment',num2str(experiment)];
                        nameExperiment = ['CTinv=',CTinv_cell{i},', CTweight=',num2str(CTweight_mat(w),'%.2f'),', R=',num2str(R_mat(r),'%.2f'),', Scale=',num2str(scale_cell{s}),', Quant.Algo=',algo_cell{a},', Ng=',num2str(Ng_mat(n))];
                        textures.List.(strExperiment) = nameExperiment;
                        fprintf(['PERFORMING EXPERIMENT %u OF %u: ''',nameExperiment,''' ... '],experiment,nExperiment)
                        [ROIonly,levels] = prepareVolume(ROIbox,maskBox,'Other',pixelW,pixelW,1,'pixelW','Matrix',algo_cell{a},Ng_mat(n)); % Pre-processing, WBPF and resampling already applied. Thus, we insert 'Other', R=1, sliceS = pixelW / Scale='pixelW', respectively
                        [GLCM] = getGLCM(ROIonly,levels); [GLCM_text] = getGLCMtextures(GLCM);
                        [GLRLM] = getGLRLM(ROIonly,levels); [GLRLM_text] = getGLRLMtextures(GLRLM);
                        [GLSZM] = getGLSZM(ROIonly,levels); [GLSZM_text] = getGLSZMtextures(GLSZM);
                        [NGTDM,countValid] = getNGTDM(ROIonly,levels); [NGTDM_text] = getNGTDMtextures(NGTDM,countValid);
                        textures.(strExperiment).Global = Global_text;
                        textures.(strExperiment).GLCM = GLCM_text;
                        textures.(strExperiment).GLRLM = GLRLM_text;
                        textures.(strExperiment).GLSZM = GLSZM_text;
                        textures.(strExperiment).NGTDM = NGTDM_text;
                        textures.(strExperiment).parameters.CTinv = CTinv_cell{i};
                        textures.(strExperiment).parameters.CTweight = CTweight_mat(w);
                        textures.(strExperiment).parameters.R = R_mat(r);
                        textures.(strExperiment).parameters.Scale = scale_cell{s};
                        textures.(strExperiment).parameters.Algo = algo_cell{a};
                        textures.(strExperiment).parameters.Ng = Ng_mat(n);
                        toc
                    end
                end
            end
        end
    end
end

end