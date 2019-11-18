function [percentAUCdecrease,varNames] = calcFeatureImportanceRFvolClin_HN(pathRF,nPerms,nameOutcome,fSetName,seed)

startpath = pwd;

cd(pathRF), load('volClinicalPARAM'), load('volumeTest'), load('testing')
if strcmp(nameOutcome,'DeathSign')
    outcome = testing.outcomes.Death;
else
    outcome = testing.outcomes.(nameOutcome);
end
nInst = numel(outcome); % Number of testing patients
perms = zeros(nInst,nPerms); rng(seed)
for p = 1:nPerms
    perms(:,p) = datasample(1:nInst,nInst,'Replace',false)';
end

results = load(['testResultsRF_',fSetName,'_',nameOutcome]); results = struct2cell(results); results = results{1};
AUC_baseline = results.AUC;
RF = load(['RF_',fSetName,'_',nameOutcome]); RF = struct2cell(RF); RF = RF{1};
indClinic = parameters.clinical.(nameOutcome);
Volume = volumeTest; tableVol = table(Volume);
tableTest = [tableVol,testing.clinical.table(:,indClinic)];
varNames = tableTest.Properties.VariableNames'; nVar = numel(varNames);
percentAUCdecrease = zeros(nVar,1);
for v = 1:nVar
    tableTemp = tableTest;
    percentAUCdecrease_temp = 0;
    for p = 1:nPerms
        tableTemp(:,v) = tableTest(perms(:,p),v);
        [prob] = predictRF(tableTemp,RF);
        [AUC,~,~,~] = calcPerformMetrics(prob,outcome,0.5);
        percentAUCdecrease_temp = percentAUCdecrease_temp + (AUC - AUC_baseline)/AUC_baseline; 
    end
    if percentAUCdecrease_temp > 0, percentAUCdecrease_temp = percentAUCdecrease_temp * -1; end
    percentAUCdecrease_temp = percentAUCdecrease_temp/nPerms;
    percentAUCdecrease(v) = percentAUCdecrease_temp;
end
[percentAUCdecrease,ind] = sort(percentAUCdecrease,'descend');
varNames = varNames(ind);

cd(startpath)
end