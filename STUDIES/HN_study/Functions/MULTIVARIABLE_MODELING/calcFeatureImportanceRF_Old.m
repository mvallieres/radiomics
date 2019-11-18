function calcFeatureImportanceRF_Old(pathRF,nPerms,seed)

startpath = pwd;

cd(pathRF), load('training'), load('testing')
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
fSetNames = fieldnames(training.textures.(nameOutcomes{1})); nFset = numel(fSetNames);

nInst = numel(testing.outcomes.(nameOutcomes{1})); % Number of testing patients
perms = zeros(nInst,nPerms);
rng(seed)
for p = 1:nPerms
    perms(:,p) = datasample(1:nInst,nInst,'Replace',false)';
end

for o = 1:nOutcomes
    outcome = testing.outcomes.(nameOutcomes{o});
    for f = 1:nFset
        results = load(['testResults_',fSetNames{f},'_',nameOutcomes{o}]); results = struct2cell(results); results = results{1};
        AUC_baseline = results.AUC;
        RF = load(['RF_',fSetNames{f},'_',nameOutcomes{o}]); RF = struct2cell(RF); RF = RF{1};
        indClinic = training.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{f});
        text = testing.textures.(nameOutcomes{o}).(fSetNames{f}); nText = size(text,2);
        tableTest = [text,testing.clinical.table(:,indClinic)];
        varNames = tableTest.Properties.VariableNames'; nVar = numel(varNames);
        percentAUCdecrease = zeros(nVar,1);
        for v = 1:nVar
            tableTemp = tableTest;
            percentAUCdecrease_temp  = 0;
            for p = 1:nPerms
                tableTemp(:,v) = tableTest(perms(:,p),v);
                [prob] = predictRF(tableTemp,RF);
                [AUC,~,~,~] = calcPerformMetrics(prob,outcome,0.5);
                percentAUCdecrease_temp = percentAUCdecrease_temp + (AUC - AUC_baseline)/AUC_baseline; 
            end
            percentAUCdecrease_temp = percentAUCdecrease_temp/nPerms;
            percentAUCdecrease(v) = percentAUCdecrease_temp;
        end
        [percentAUCdecrease,ind] = sort(percentAUCdecrease,'descend');
        variableImportance.(nameOutcomes{o}).(fSetNames{f}).percentAUCdecrease = percentAUCdecrease;
        varNames = varNames(ind);
        variableImportance.(nameOutcomes{o}).(fSetNames{f}).varNames = varNames;
    end
end
save('testingVariableImportance','variableImportance')

cd(startpath)
end