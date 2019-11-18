function calcFeatureImportanceRFvolClin_batchHN(pathRF,nPerms,matlabPATH,seed)

startpath = pwd;

time = 60;
cd(pathRF), load('training')
mkdir('batchLog_RFpermVolClin'), cd('batchLog_RFpermVolClin'), pathBatch = pwd;
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
% fSetNames = fieldnames(training.textures.(nameOutcomes{1})); nFset = numel(fSetNames);
fSetNames = {'volClinical'}; nFset = numel(fSetNames);


% PRODUCE BATCH COMPUTATIONS
batch = 0;
save('workspace','pathRF','nPerms','nameOutcomes','fSetNames','seed')
for o = 1:nOutcomes
    for f = 1:nFset
        batch = batch + 1;
        nameScript = ['batch',num2str(batch),'_script.m'];
        fid = fopen(nameScript,'w');
        fprintf(fid,'load(''workspace'')\n');
        fprintf(fid,['[percentAUCdecrease,varNames] = calcFeatureImportanceRFvolClin_HN(pathRF,nPerms,nameOutcomes{',num2str(o),'},fSetNames{',num2str(f),'},seed);\n']);
        fprintf(fid,['save(''result_batch',num2str(batch),''',''percentAUCdecrease'',''varNames'')\n']);
        fprintf(fid,['system(''touch batch',num2str(batch),'_end'');\n']);
        fprintf(fid,'clear all');
        fclose(fid);
        system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
    end
end

% WAITING LOOP
waitBatch(pathBatch,time,nOutcomes*nFset)
delete('workspace.mat')


% GROUP RESULTS IN THE 'training' STRUCT
cd(pathRF)
if ~exist('testingVariableImportance.mat')
    variableImportance = struct;
else
    load('testingVariableImportance')
end
cd(pathBatch)
batch = 0;
for o = 1:nOutcomes
    for f = 1:nFset
        batch = batch + 1;
        load(['result_batch',num2str(batch)]) % Variables 'percentAUCdecrease' and 'varNames' gets out of there
        variableImportance.(nameOutcomes{o}).(fSetNames{f}).percentAUCdecrease = percentAUCdecrease;
        variableImportance.(nameOutcomes{o}).(fSetNames{f}).varNames = varNames;
        delete(['result_batch',num2str(batch),'.mat'])
        clear percentAUCdecrease varNames
    end
end
cd(pathRF), save('testingVariableImportance','variableImportance')

cd(startpath)
end