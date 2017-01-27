function performBatchRFtextClin_LGG(pathRF,nTrees,seed,matlabPATH)
% -------------------------------------------------------------------------
% function performBatchRFtextClin_LGG(pathRF,nTrees,seed,matlabPATH)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function performs prediction estimation using out-of-bag estimates
% for random forests constructed using (TEXTURE + CLINICAL) features.
% -------------------------------------------------------------------------
% INPUTS:
% 1. pathRF: Full path to the RANDOM_FORESTS results folder.
%            --> Ex: '/myProject/WORKSPACE/RANDOM_FORESTS'
% 2. nTrees: Number of decision-trees in the random forests.
%            --> Ex: 500
% 3. seed: Numerical number to use as seed for random-forest construction.
%          --> Ex: 54288
% 6. matlabPATH: Full path to the MATLAB executable on the system.
%                --> 'matlab' if a symbolic link to the matlab executable
%                     was previously created. 
% -------------------------------------------------------------------------
% OUTPUTS: Final out-of-bag prediction estimation results saved in a table
% in 'pathRF'.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2017
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015-2017  Martin Vallieres
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

cd(pathRF), load('dataStruct')
nameOutcomes = fieldnames(dataStruct); nOutcomes = numel(nameOutcomes);

time = 30; testCost = 0.25:0.25:5;
mkdir('batchLog_RFtextClin'), cd('batchLog_RFtextClin'), pathBatch = pwd;


% PRODUCE BATCH COMPUTATIONS
save('workspace','nTrees','dataStruct','nameOutcomes','seed','testCost'), batch = 0;
for o = 1:nOutcomes
    batch = batch + 1;
    nameScript = ['batch',num2str(batch),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['nameOutcome = nameOutcomes{',num2str(o),'};\n']);
    fprintf(fid,['nVas = dataStruct.(nameOutcome).nVas;\n']);
    fprintf(fid,['nText = dataStruct.(nameOutcome).nText;\n']);
    fprintf(fid,['X = dataStruct.(nameOutcome).data(:,(nVas+1):end);\n']);
    fprintf(fid,['Y = dataStruct.(nameOutcome).outcome;\n']);
    fprintf(fid,['cat = dataStruct.(nameOutcome).categories((nVas+1):end);\n']);
    fprintf(fid,['[auc,sensitivity,specificity,accuracy] = predictionEstimationRF_LGG(X,Y,logical(cat),nTrees,seed,testCost);\n']);
    fprintf(fid,['save(''result_batch',num2str(batch),''',''auc'',''sensitivity'',''specificity'',''accuracy'')\n']);
    fprintf(fid,['system(''touch batch',num2str(batch),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nOutcomes)
delete('workspace.mat')

% GROUP RESULTS
cd(pathBatch)
batch = 0; AUC = zeros(nOutcomes,1); Sensitivity = zeros(nOutcomes,1); Specificity = zeros(nOutcomes,1); Accuracy = zeros(nOutcomes,1);
for o = 1:nOutcomes
    batch = batch + 1;
    load(['result_batch',num2str(batch)]) % Variable 'auc', 'sensitivity', 'specificity', 'accuracy' gets out of here.
    AUC(o) = auc; Sensitivity(o) = sensitivity; Specificity(o) = specificity; Accuracy(o) = accuracy;
    delete(['result_batch',num2str(batch),'.mat'])
    clear auc sensitivity specificity accuracy
end
tab = table(AUC,Sensitivity,Specificity,Accuracy,'RowNames',nameOutcomes);
cd(pathRF), save('tableSummary_TextClin','tab')

cd(startpath)
end