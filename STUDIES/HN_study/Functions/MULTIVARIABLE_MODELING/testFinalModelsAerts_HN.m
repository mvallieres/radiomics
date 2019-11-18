function testFinalModelsAerts_HN(pathExperiment,fSetNames)
% - paramCells: One cell of parameters for each fSet
% - nameCells: One cell of names for each fSet

startpath = pwd;

cd(pathExperiment), load('training'), load('testing')
cd('FINAL_MODELS'), pathFinalModels = pwd;
nameOutcomes = fieldnames(training); nOutcomes = numel(nameOutcomes); 
nFset = numel(fSetNames);

for o = 1:nOutcomes
    for f = 1:nFset
        results = []; results = struct;
        cd(fullfile(pathFinalModels,nameOutcomes{o},[fSetNames{f},'sign']))
        load('finalModel'), load('coeff'), load('response'), order = 4;
        results.model.Name = finalModel.Name; results.model.coeff = coeff; results.model.order = order;
        results.trainData.data = finalModel.Data; results.trainData.response = response; results.trainData.outcome = training.(nameOutcomes{o}).outcome;

        % TESTING MODEL
        outcome = testing.(nameOutcomes{o}).outcome;
        resp = zeros(numel(outcome),1);
        data = zeros(numel(outcome),order);
        for n = 1:order
            data(:,n) = testing.(nameOutcomes{o}).sign.(fSetNames{f})(:,n);
            resp = resp + coeff(n)*data(:,n);
        end
        resp = resp + coeff(end);
        results.testData.data = data; results.testData.response = resp; results.testData.outcome = outcome;
        [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(resp,outcome,0);
        results.AUC = AUC;
        results.Sensitivity = sensitivity;
        results.Specificity = specificity;
        results.Accuracy = accuracy;
        cd(pathExperiment), save(['testResultsLR_',[fSetNames{f},'sign'],'_',nameOutcomes{o}],'results')
    end
end

cd(startpath)
end