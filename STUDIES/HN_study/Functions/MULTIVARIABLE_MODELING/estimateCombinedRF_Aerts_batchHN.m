function estimateCombinedRF_Aerts_batchHN(pathRF,nSplit,nBoot,testSplit,testCost,seed,matlabPATH)

startpath = pwd;

% INITIALIZATIONS
cd(pathRF)
load('clinicalPARAM') % parameters gets out of there
training = load('training'); training = struct2cell(training); training = training{1};
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
fSetNames = {'PETsign','CTsign','CTorigsign'}; nFset = numel(fSetNames); % The three sets tested for the radiomic signature
comb = cell(1,nOutcomes);
for o = 1:nOutcomes
    comb{o} = {parameters.clinical.(nameOutcomes{o})};
end
time = 60;
mkdir('batchLog_RFsign'), cd('batchLog_RFsign'), pathBatch = pwd;


% PRODUCE BATCH COMPUTATIONS FOR RF ESTIMATION
batch = 0;
save('workspace','pathRF','nameOutcomes','fSetNames','comb','testCost','nSplit','testSplit','nBoot','seed'), pause(2);
for o = 1:nOutcomes
    for f = 1:nFset
        batch = batch + 1;
        nameScript = ['batch',num2str(batch),'_script.m'];
        fid = fopen(nameScript,'w');
        fprintf(fid,'load(''workspace'')\n');
        fprintf(fid,['[indClinic,bestCost] = estimateCombinedRF(pathRF,nameOutcomes{',num2str(o),'},fSetNames{',num2str(f),'},comb{',num2str(o),'},testCost,nSplit,testSplit,nBoot,seed);\n']);
        fprintf(fid,['save(''result_batch',num2str(batch),''',''indClinic'',''bestCost'')\n']);
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
cd(pathBatch)
batch = 0;
for o = 1:nOutcomes
    for f = 1:nFset
        batch = batch + 1;
        load(['result_batch',num2str(batch)]) % Variables 'indClinic' and 'bestCost' gets out of there
        training.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{f}) = indClinic;
        training.cost.(nameOutcomes{o}).(fSetNames{f}) = bestCost;
        delete(['result_batch',num2str(batch),'.mat'])
        clear indClinic bestCost
    end
end
cd(pathRF), save('training','training')

cd(startpath)
end