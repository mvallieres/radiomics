function testFinalModelsTimeAerts_HN(pathExperiment,fSetNames)
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
        load('finalModel'), load('coeff'), load('medianHR'), load('response'), order = 4;
        results.model.Name = finalModel.Name; results.model.coeff = coeff; results.model.medianHR = medianHR; results.model.order = order;
        results.trainData.data = finalModel.Data; results.trainData.response = response; results.trainData.outcome = training.Death.outcome; results.trainData.timeToEvent = training.Death.timeToEvent; results.trainData.censoring = 1 - results.trainData.outcome;
        
        % TESTING THE MODEL
        time = testing.Death.timeToEvent;
        censoring = 1 - testing.Death.outcome;
        resp = zeros(numel(time),1);
        data = zeros(numel(time),order);
        for n = 1:order
            data(:,n) = testing.(nameOutcomes{o}).sign.(fSetNames{f})(:,n);
            resp = resp + coeff(n)*data(:,n);
        end
        results.testData.data = data; results.testData.response = resp; results.testData.outcome = 1 - censoring; results.testData.timeToEvent = time; results.testData.censoring = censoring;
        CI = calcCI(resp,time,censoring);
        results.CI = CI;
        cd(pathExperiment), save(['testResultsCR_',[fSetNames{f},'sign'],'_',nameOutcomes{o}],'results')
    end
end

cd(startpath)
end