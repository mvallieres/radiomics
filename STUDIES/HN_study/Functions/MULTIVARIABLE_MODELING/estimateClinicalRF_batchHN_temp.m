function estimateClinicalRF_batchHN_temp(pathRF,nSplit,nBoot,matlabPATH,seed)

startpath = pwd;

% INITIALIZATIONS
cd(pathRF)
training = load('training'); training = struct2cell(training); training = training{1};
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
%fSetNames = fieldnames(training.textures.(nameOutcomes{1})); nFset = numel(fSetNames);
fSetNames = {'PET','CT','CT'}; nFset = 1;% Here, only 'PET' is tested for Loco, 'CT for 'Distant' and 'CT' for 'Death' 'f' was replacef by 'o' in three lines)
clinicNames = training.clinical.table.Properties.VariableNames';
testSplit = 1/3; % Proportion of test cases in the splits
comb =  {{[1,2,3],...      % Loco       : Age, Subtype, T_Stage
        [1,2,4],...        % Loco       : Age, Subtype, T_Stage, N_Stage
        [1,2,5],...        % Loco       : Age, Subtype, TNM_Stage
        [1,2,3,4]},...     % Loco       : Age, Subtype, T_Stage, N_Stage
        {[1,2,3],...       % Distant       : Age, Subtype, T_Stage
        [1,2,4],...        % Distant       : Age, Subtype, N_Stage
        [1,2,5],...        % Distant       : Age, Subtype, TNM_Stage
        [1,2,3,4]},...     % Distant       : Age, Subtype, T_Stage, N_Stage
        {[1,2,3],...       % Death       : Age, Subtype, T_Stage
        [1,2,4],...        % Death       : Age, Subtype, N_Stage
        [1,2,5],...        % Death       : Age, Subtype, TNM_Stage
        [1,2,3,4]}};       % Death       : Age, Subtype, T_Stage, N_Stage
    
nComb = [4,4,4];
testCost = 0.4:0.1:2.5; % Emphasis factor on positive instances during random forest training

time = 60;
mkdir('batchLog_RF'), cd('batchLog_RF'), pathBatch = pwd;


% PRODUCE BATCH COMPUTATIONS (for each outcome/fSet combination, 3 in total)
save('workspace','pathRF','nameOutcomes','fSetNames','comb','testCost','nSplit','testSplit','nBoot','seed'), pause(2);
for i = 1:nOutcomes
    nameScript = ['batch',num2str(i),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['[indClinic,bestCost] = estimateClinicalRF(pathRF,nameOutcomes{',num2str(i),'},fSetNames{',num2str(i),'},comb{',num2str(i),'},testCost,nSplit,testSplit,nBoot,seed);\n']);
    fprintf(fid,['save(''result_batch',num2str(i),''',''indClinic'',''bestCost'')\n']);
    fprintf(fid,['system(''touch batch',num2str(i),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nOutcomes)
delete('workspace.mat')


% GROUP RESULTS IN THE 'training' STRUCT
cd(pathBatch)
for o = 1:nOutcomes
    load(['result_batch',num2str(o)]) % Variables 'indClinic' and 'bestCost' gets out of there
    training.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{o}) = indClinic;
    training.cost.(nameOutcomes{o}).(fSetNames{o}) = bestCost;
    delete(['result_batch',num2str(o),'.mat'])
    clear indClinic bestCost
end
cd(pathRF), save('training','training')

cd(startpath)
end