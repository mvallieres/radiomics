function computeVolumeResults_Regression(pathLR,pathCR,imbalance,seed,pathFig)

startpath = pwd;

% LOGISTIC REGRESSION RESULTS
cd(pathLR), load('training'), load('testing')
nameOutcomes = fieldnames(training); nOutcomes = numel(nameOutcomes);
for o = 1:nOutcomes
    Xtrain = training.(nameOutcomes{o}).nonText.Volume.Data;
    Ytrain = training.(nameOutcomes{o}).outcome;
    Xtest = testing.(nameOutcomes{o}).nonText.Volume.Data;
    Ytest = testing.(nameOutcomes{o}).outcome;
    [coeff,respTrain,~] = computeModelCoefficients_HN(Xtrain,Ytrain,imbalance,seed);
    [respTest] = responseLR(Xtest,coeff);
    [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(respTest,Ytest,0);
    results = struct;
    results.model.Feature1.isTexture = false; results.model.Feature1.nonText = 'Volume'; results.model.coeff = coeff; results.model.order = 1;
    results.trainData.data = Xtrain; results.trainData.response = respTrain; results.trainData.outcome = Ytrain;
    results.testData.data = Xtest; results.testData.response = respTest; results.testData.outcome = Ytest;
    results.AUC = AUC; results.Sensitivity = sensitivity; results.Specificity = specificity; results.Accuracy = accuracy;
    save(['testResultsLR_VOL_',nameOutcomes{o}],'results'), clear results
end


% COX REGRESSION RESULTS + KAPLAN-MEIER CURVES
cd(pathCR), load('training'), load('testing')
nameOutcomes = fieldnames(training); nOutcomes = numel(nameOutcomes);
for o = 1:nOutcomes
    Xtrain = training.(nameOutcomes{o}).nonText.Volume.Data;
    Ytrain = training.(nameOutcomes{o}).timeToEvent;
    censTrain = 1 - training.(nameOutcomes{o}).outcome;
    Xtest = testing.(nameOutcomes{o}).nonText.Volume.Data;
    Ytest = testing.(nameOutcomes{o}).timeToEvent;
    censTest = 1 - testing.(nameOutcomes{o}).outcome;
    [coeff,respTrain,~,medianHR] = computeModelCoefficientsTime_HN(Xtrain,Ytrain,censTrain,seed);
    [respTest] = responseCox(Xtest,coeff);
    CI = calcCI(respTest,Ytest,censTest);
    results = struct;
    results.model.Feature1.isTexture = false; results.model.Feature1.nonText = 'Volume'; results.model.coeff = coeff; results.model.medianHR = medianHR; results.model.order = 1;
    results.trainData.data = Xtrain; results.trainData.response = respTrain; results.trainData.outcome = 1-censTrain; results.trainData.timeToEvent = Ytrain; results.trainData.censoring = censTrain;
    results.testData.data = Xtest; results.testData.response = respTest; results.testData.outcome = 1 - censTest; results.testData.timeToEvent = Ytest; results.testData.censoring = censTest;
    results.CI = CI;
    save(['testResultsCR_VOL_',nameOutcomes{o}],'results'), clear results
    displayKaplanMeierCR_HN(pathCR,{nameOutcomes{o}},{'VOL'},pathFig)
end

cd(startpath)
end