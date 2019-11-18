function computeClinicalResults_RF(pathRF,nBoot,nSplit,testSplit,testCost,seeds,matlabPATH)

startpath = pwd;

cd(pathRF), load('training'), load('testing')
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
comb =  {{[1,2,3],...      % Loco        : Age, Subtype, T_Stage
        [1,2,4],...        % Loco        : Age, Subtype, N_Stage
        [1,2,5],...        % Loco        : Age, Subtype, TNM_Stage
        [1,2,3,4]},...     % Loco        : Age, Subtype, T_Stage, N_Stage
        {[1,2,3],...       % Distant     : Age, Subtype, T_Stage
        [1,2,4],...        % Distant     : Age, Subtype, N_Stage
        [1,2,5],...        % Distant     : Age, Subtype, TNM_Stage
        [1,2,3,4]},...     % Distant     : Age, Subtype, T_Stage, N_Stage
        {[1,2,3],...       % Death       : Age, Subtype, T_Stage
        [1,2,4],...        % Death       : Age, Subtype, N_Stage
        [1,2,5],...        % Death       : Age, Subtype, TNM_Stage
        [1,2,3,4]}};       % Death       : Age, Subtype, T_Stage, N_Stage
mkdir('batchLog_RFclinical'), cd('batchLog_RFclinical'), pathBatch = pwd;
time = 60; % Time in the wainting loop

% PRODUCE BATCH COMPUTATIONS TO FIND THE BEST CLINICAL VARIABLES
seed = seeds(1); batch = 0;
save('workspace','pathRF','nameOutcomes','comb','testCost','nSplit','testSplit','nBoot','seed'), pause(2);
for o = 1:nOutcomes
    batch = batch + 1;
    nameScript = ['batch',num2str(batch),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['[indClinic,bestCost] = estimateClinicalRF(pathRF,nameOutcomes{',num2str(o),'},comb{',num2str(o),'},testCost,nSplit,testSplit,nBoot,seed);\n']);
    fprintf(fid,['save(''result_batch',num2str(batch),''',''indClinic'',''bestCost'')\n']);
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
batch = 0; parameters = struct;
for o = 1:nOutcomes
    batch = batch + 1;
    load(['result_batch',num2str(batch)]) % Variables 'indClinic' and 'bestCost' gets out of there
    parameters.clinical.(nameOutcomes{o}) = indClinic;
    parameters.cost.(nameOutcomes{o}) = bestCost;
    delete(['result_batch',num2str(batch),'.mat'])
    clear indClinic bestCost
end
cd(pathRF), save('clinicalPARAM','parameters')


% COMPUTE THE RANDOM FORESTS
for o = 1:nOutcomes
    if strcmp(nameOutcomes{o},'DeathSign')
        outcome = training.outcomes.Death;
    else
        outcome = training.outcomes.(nameOutcomes{o});
    end
    indClinic = parameters.clinical.(nameOutcomes{o});
    cost = parameters.cost.(nameOutcomes{o});
    tableTrain = training.clinical.table(:,indClinic);
    cat = logical(training.clinical.categories(indClinic));
    rng(seeds(2)), [RF] = trainRF_table(tableTrain,outcome,cat,nBoot,cost);
    RF = compact(RF); % Compact version
    save(['RF_clinic_',nameOutcomes{o}],'RF')
end


% TEST THE RANDOM FORESTS
cd(pathRF), load('clinicalPARAM')
for o = 1:nOutcomes
    if strcmp(nameOutcomes{o},'DeathSign')
        outcome = testing.outcomes.Death;
        time = testing.timeToEvents.Death;
    else
        outcome = testing.outcomes.(nameOutcomes{o});
        time = testing.timeToEvents.(nameOutcomes{o});
    end
    censoring = 1 - outcome;
    results = struct;
    RF = load(['RF_clinic_',nameOutcomes{o}]); RF = struct2cell(RF); RF = RF{1};
    indClinic = parameters.clinical.(nameOutcomes{o});
    tableTest = testing.clinical.table(:,indClinic);
    [prob] = predictRF(tableTest,RF);
    results.probResponse = prob;
    [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(prob,outcome,0.5);
    CI = calcCI(prob,time,censoring);
    results.AUC = AUC; results.Sensitivity = sensitivity; results.Specificity = specificity; results.Accuracy = accuracy;
    results.CI = CI;
    save(['testResultsRF_clinic_',nameOutcomes{o}],'results'), clear results
end

cd(startpath)
end