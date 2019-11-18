function computeVolumeResults_RF(pathRF,nBoot,nSplit,testSplit,testCost,seeds,matlabPATH)

startpath = pwd;

cd(pathRF), load('training'), load('testing'), load('clinicalPARAM')
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
comb = cell(1,nOutcomes);
for o = 1:nOutcomes
    comb{o} = {parameters.clinical.(nameOutcomes{o})};
end
time = 60;
mkdir('batchLog_RFvolClin'), cd('batchLog_RFvolClin'), pathBatch = pwd;


% PRODUCE BATCH COMPUTATIONS FOR RF ESTIMATION
seed = seeds(1); batch = 0;
save('workspace','pathRF','nameOutcomes','comb','testCost','nSplit','testSplit','nBoot','seed'), pause(2);
for o = 1:nOutcomes
    batch = batch + 1;
    nameScript = ['batch',num2str(batch),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['[indClinic,bestCost] = estimateVolClinicalRF(pathRF,nameOutcomes{',num2str(o),'},comb{',num2str(o),'},testCost,nSplit,testSplit,nBoot,seed);\n']);
    fprintf(fid,['save(''result_batch',num2str(batch),''',''indClinic'',''bestCost'')\n']);
    fprintf(fid,['system(''touch batch',num2str(batch),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP FOR RF ESTIMATION
waitBatch(pathBatch,time,nOutcomes)
delete('workspace.mat')

% GROUP RESULTS FOR RF ESTIMATION
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
cd(pathRF), save('volClinicalPARAM','parameters')


% COMPUTE THE FINAL RANDOM FORESTS
load('volumeTrain'), Volume = volumeTrain; tableVol = table(Volume);
for o = 1:nOutcomes
    if strcmp(nameOutcomes{o},'DeathSign')
        outcome = training.outcomes.Death;
    else
        outcome = training.outcomes.(nameOutcomes{o});
    end
    indClinic = parameters.clinical.(nameOutcomes{o});
    cost = parameters.cost.(nameOutcomes{o});
    tableTrain = [tableVol,training.clinical.table(:,indClinic)];
    cat = [false,logical(training.clinical.categories(indClinic))];
    rng(seeds(2)), [RF] = trainRF_table(tableTrain,outcome,cat,nBoot,cost);
    RF = compact(RF); % Compact version
    save(['RF_VOLclinic_',nameOutcomes{o}],'RF')
end


% TEST THE FINAL RANDOM FORESTS
load('volumeTest'), Volume = volumeTest; tableVol = table(Volume);
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
    RF = load(['RF_VOLclinic_',nameOutcomes{o}]); RF = struct2cell(RF); RF = RF{1};
    indClinic = parameters.clinical.(nameOutcomes{o});
    tableTest = [tableVol,testing.clinical.table(:,indClinic)];
    [prob] = predictRF(tableTest,RF);
    results.probResponse = prob;
    [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(prob,outcome,0.5);
    CI = calcCI(prob,time,censoring);
    results.AUC = AUC; results.Sensitivity = sensitivity; results.Specificity = specificity; results.Accuracy = accuracy;
    results.CI = CI;
    save(['testResultsRF_VOLclinic_',nameOutcomes{o}],'results'), clear results
end

cd(startpath)
end