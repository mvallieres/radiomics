function estimateClinicalRF_batchHN(pathRF,nSplit,nBoot,matlabPATH,seed)

startpath = pwd;

% INITIALIZATIONS
cd(pathRF)
training = load('training'); training = struct2cell(training); training = training{1};
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
clinicNames = training.clinical.table.Properties.VariableNames';
testSplit = 1/3; % Proportion of test cases in the splits
comb =  {{[1,2,3],...      % Loco        : Age, Subtype, T_Stage
        [1,2,5],...        % Loco        : Age, Subtype, TNM_Stage
        [1,2,3,4]},...     % Loco        : Age, Subtype, T_Stage, N_Stage
        {[1,2,4],...       % Distant     : Age, Subtype, N_Stage
        [1,2,5],...        % Distant     : Age, Subtype, TNM_Stage
        [1,2,3,4]},...     % Distant     : Age, Subtype, T_Stage, N_Stage
        {[1,2,3],...       % Death       : Age, Subtype, T_Stage
        [1,2,4],...        % Death       : Age, Subtype, N_Stage
        [1,2,5],...        % Death       : Age, Subtype, TNM_Stage
        [1,2,3,4]}};       % Death       : Age, Subtype, T_Stage, N_Stage
testCost = 0.4:0.1:2.5; % Emphasis factor on positive instances during random forest training

time = 60;
mkdir('batchLog_RF'), cd('batchLog_RF'), pathBatch = pwd;
fSetNames = {'PET','CT','PETCT'}; nFset = numel(fSetNames); % TESTING PET radiomics models, CT radiomics models, and PET radiomics models + CT radiomics models, as organized in organizeRFexperiments.m


% PRODUCE BATCH COMPUTATIONS (for each outcome/fSet combination, 3 in total)
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