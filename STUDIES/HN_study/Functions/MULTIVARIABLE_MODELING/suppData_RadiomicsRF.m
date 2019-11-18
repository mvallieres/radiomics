function suppData_RadiomicsRF(pathSUPP,nBoot,nSplit,testSplit,testCost,seeds,matlabPATH,pathFig)

startpath = pwd;

cd(pathSUPP), load('training'), load('testing')
nameOutcomes = {'Locoregional','Distant','Death'}; nOutcomes = numel(nameOutcomes);
fSetNames = {'PETCT','CT','PET'};
nExp = numel(nameOutcomes);
time = 60;
mkdir('batchLog_RFradiomics'), cd('batchLog_RFradiomics'), pathBatch = pwd;

% PRODUCE BATCH COMPUTATIONS (for each outcome/fSet combination, 3 in total)
seed = seeds(1); batch = 0;
save('workspace','pathSUPP','nameOutcomes','fSetNames','testCost','nSplit','testSplit','nBoot','seed'), pause(2);
for o = 1:nOutcomes
    batch = batch + 1;
    nameScript = ['batch',num2str(batch),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['[bestCost] = estimateRadiomicsRF(pathSUPP,nameOutcomes{',num2str(o),'},fSetNames{',num2str(o),'},testCost,nSplit,testSplit,nBoot,seed);\n']);
    fprintf(fid,['save(''result_batch',num2str(batch),''',''bestCost'')\n']);
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
    load(['result_batch',num2str(batch)]) % Variable 'bestCost' gets out of there
    parameters.cost.(nameOutcomes{o}).(fSetNames{o}) = bestCost;
    delete(['result_batch',num2str(batch),'.mat'])
    clear indClinic bestCost
end
cd(pathSUPP), save('radiomicsPARAM','parameters')


% COMPUTE THE RANDOM FORESTS
for o = 1:nOutcomes
    if strcmp(nameOutcomes{o},'DeathSign')
        outcome = training.outcomes.Death;
    else
        outcome = training.outcomes.(nameOutcomes{o});
    end
    text = training.textures.(nameOutcomes{o}).(fSetNames{o}); nText = size(text,2);
    cost = parameters.cost.(nameOutcomes{o}).(fSetNames{o});
    tableTrain = text;
    cat = logical([zeros(1,nText)]);
    rng(seeds(2)), [RF] = trainRF_table(tableTrain,outcome,cat,nBoot,cost);
    RF = compact(RF); % Compact version
    save(['RFradiomics_',fSetNames{o},'_',nameOutcomes{o}],'RF')
end


% TEST THE RANDOM FORESTS
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
    RF = load(['RFradiomics_',fSetNames{o},'_',nameOutcomes{o}]); RF = struct2cell(RF); RF = RF{1};
    text = testing.textures.(nameOutcomes{o}).(fSetNames{o});
    tableTest = text;
    [prob] = predictRF(tableTest,RF);
    results.probResponse = prob;
    [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(prob,outcome,0.5);
    CI = calcCI(prob,time,censoring);
    results.AUC = AUC; results.Sensitivity = sensitivity; results.Specificity = specificity; results.Accuracy = accuracy;
    results.CI = CI;
    save(['testResultsRFradiomics_',fSetNames{o},'_',nameOutcomes{o}],'results'), clear results
end


% DISPLAY KAPLAN-MEIER CURVES
displayKaplanMeierRF_HNsupp(pathSUPP,{'Locoregional'},['radiomics_',fSetNames{1}],pathFig)
displayKaplanMeierRF_HNsupp(pathSUPP,{'Distant'},['radiomics_',fSetNames{2}],pathFig)
displayKaplanMeierRF_HNsupp(pathSUPP,{'Death'},['radiomics_',fSetNames{3}],pathFig)

cd(startpath)
end