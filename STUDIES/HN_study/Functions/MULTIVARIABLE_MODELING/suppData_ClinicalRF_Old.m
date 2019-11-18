function suppData_ClinicalRF_Old(pathSUPP,nBoot,seed)

startpath = pwd;

cd(pathSUPP), load('training'), load('testing')

nameOutcomes = {'Locoregional','Distant','Death','Death','DeathSign'};
fSetNames = {'PET','CT','PETCT','CT','CT'};
nExp = numel(nameOutcomes);

% COMPUTE THE RANDOM FORESTS
for o = 1:nExp
    if strcmp(nameOutcomes{o},'DeathSign')
        outcome = training.outcomes.Death;
    else
        outcome = training.outcomes.(nameOutcomes{o});
    end
    indClinic = training.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{o});
    cost = training.cost.(nameOutcomes{o}).(fSetNames{o});
    tableTrain = training.clinical.table(:,indClinic);
    cat = logical(training.clinical.categories(indClinic));
    rng(seed), [RF] = trainRF_table(tableTrain,outcome,cat,nBoot,cost);
    RF = compact(RF); % Compact version
    save(['RFclinical_',fSetNames{o},'_',nameOutcomes{o}],'RF')
end


% TEST THE RANDOM FORESTS
for o = 1:nExp
    if strcmp(nameOutcomes{o},'DeathSign')
        outcome = testing.outcomes.Death;
        time = testing.timeToEvents.Death;
    else
        outcome = testing.outcomes.(nameOutcomes{o});
        time = testing.timeToEvents.(nameOutcomes{o});
    end
    censoring = 1 - outcome;
    results = struct;
    RF = load(['RFclinical_',fSetNames{o},'_',nameOutcomes{o}]); RF = struct2cell(RF); RF = RF{1};
    indClinic = training.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{o});
    tableTest = testing.clinical.table(:,indClinic);
    [prob] = predictRF(tableTest,RF);
    results.probResponse = prob;
    [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(prob,outcome,0.5);
    CI = calcCI(prob,time,censoring);
    results.AUC = AUC; results.Sensitivity = sensitivity; results.Specificity = specificity; results.Accuracy = accuracy;
    results.CI = CI;
    save(['testResultsRFclinical_',(fSetNames{o}),'_',nameOutcomes{o}],'results'), clear results
end

cd(startpath)
end