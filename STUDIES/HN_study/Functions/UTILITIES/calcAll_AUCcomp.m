function calcAll_AUCcomp(pathLR,pathRF,nameOutcomes,fSetNames)
startpath = pwd;

nOutcomes = numel(nameOutcomes);
nFset = numel(fSetNames);

% COMPUTATION
aucComp = struct;
for o = 1:nOutcomes
    for f = 1:nFset
        cd(pathLR)
        load(['testResultsLR_',fSetNames{f},'_',nameOutcomes{o}]) % results gets out of there
        X1 = logit(results.testData.response);
        Y = results.testData.outcome;
        cd(pathRF)
        load(['testResultsRF_',fSetNames{f},'clinic_',nameOutcomes{o}]) % results gets out of there
        X2 = results.probResponse;
        [p,CI] = testAUCincrease(X1,X2,Y);
        aucComp.(nameOutcomes{o}).(fSetNames{f}).p = p;
        aucComp.(nameOutcomes{o}).(fSetNames{f}).CI = CI;
    end
end

cd(pathRF), save('aucComp','aucComp')

cd(startpath)
end