function [textures] = calcPatientSepText_STS(sData,roiNumb,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% function [textures] = calcPatientSepText_STS(sData,roiNumb,R_mat,scale_cell,algo_cell,Ng_mat)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes SEPARATE textures for a single patient, for all 
% different combinations of the following texture extraction parameters:
% - WPBF ratio 'R': Ratio specifying the ratio of weight given to band-pass 
%                   coefficients to the weigth given to the rest of 
%                   coefficients (HHH and LLL) in the wavelet band-pass 
%                   filtering process.
% - Scale: Resolution at which 'volume' is isotropically resampled.
% - Quant. Algo.: Quantization algorithm used on 'volume'.
% - Ng: Number of gray-levels on the quantization process. 
%
% Different extraction parameters are passed as arrays or cells in the
% function in order to test all possible combinations. This function is 
% used for SEPARATE scans specifically. See Ref. [1] and 'prepareVolume.m' 
% for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% INPUTS:
% - sData: File organizing the data to analyze.
% - roiNumb: Contour Number to use in sData{2}.scan.contour
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
% R_mat = [1/2,2/3,1,3/2,2];
% scale_cell = {'pixelW',1,2,3,4,5};
% algo_cell = {'Equal','Lloyd'};
% Ng_mat = [8,16,32,64];
% [textures] = calcPatientSepText_STS(sData,1,R_mat,scale_cell,algo_cell,Ng_mat);
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
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
textures.Parameters.R = R_mat;
textures.Parameters.Scale = scale_cell;
textures.Parameters.Algo = algo_cell;
textures.Parameters.Ng = Ng_mat;
nExperiment = numel(R_mat)*numel(scale_cell)*numel(algo_cell)*numel(Ng_mat);
volume = getROIbox(sData,roiNumb); mask = sData{2}.scan.contour(roiNumb).boxMask;
pixelW = sData{2}.scan.pixelW; sliceS = sData{2}.scan.sliceS; scanType = sData{2}.type;

experiment = 0;
fprintf('\n')
for r = 1:numel(R_mat)
    for s = 1:numel(scale_cell)
        [ROIonly,~,ROIbox,maskBox] = prepareVolume(volume,mask,scanType,pixelW,sliceS,R_mat(r),scale_cell{s},'Global');
        [Global_text] = getGlobalTextures(ROIonly,100);
        for a = 1:numel(algo_cell)
            for n = 1:numel(Ng_mat)
                tic
                experiment = experiment + 1;
                strExperiment = ['Experiment',num2str(experiment)];
                nameExperiment = ['R=',num2str(R_mat(r),'%.2f'),', Scale=',num2str(scale_cell{s}),', Quant.Algo=',algo_cell{a},', Ng=',num2str(Ng_mat(n))];
                textures.List.(strExperiment) = nameExperiment;
                fprintf(['PERFORMING EXPERIMENT %u OF %u: ''',nameExperiment,''' ... '],experiment,nExperiment)
                [ROIonly,levels] = prepareVolume(ROIbox,maskBox,'Other',pixelW,pixelW,1,'pixelW','Matrix',algo_cell{a},Ng_mat(n)); % Pre-processing, WBPF and resampling already applied. Thus, we insert 'Other', R=1, sliceS = pixelW / Scale='pixelW', respectively
                [GLCM] = getGLCM(ROIonly,levels); [GLCM_text] = getGLCMtextures(GLCM);
                [GLRLM] = getGLRLM(ROIonly,levels); [GLRLM_text] = getGLRLMtextures_STS(GLRLM);
                [GLSZM] = getGLSZM(ROIonly,levels); [GLSZM_text] = getGLSZMtextures_STS(GLSZM);
                [NGTDM,countValid] = getNGTDM(ROIonly,levels); [NGTDM_text] = getNGTDMtextures(NGTDM,countValid);
                textures.(strExperiment).Global = Global_text;
                textures.(strExperiment).GLCM = GLCM_text;
                textures.(strExperiment).GLRLM = GLRLM_text;
                textures.(strExperiment).GLSZM = GLSZM_text;
                textures.(strExperiment).NGTDM = NGTDM_text;
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