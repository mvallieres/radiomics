function calcCorrRFprobTime(pathRF,nameOutcomes,nameSets)

startpath = pwd;
nOutcomes = numel(nameOutcomes);


cd(pathRF), load('testing')
correlations = struct;
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o}; fSet = nameSets{o};
    load(['testResults_',fSet,'_',nameOutcome]) % results gets out of there
    probData = results.probResponse;
    time = testing.timeToEvents.(nameOutcome);
    outcome = testing.outcomes.(nameOutcome);
%     probPos = probData(outcome == 1);
%     timePos = time(outcome == 1);
%     [rs,p] = corr(probPos,timePos,'type','Spearman');
[rs,p] = corr(probData,time,'type','Spearman');
    correlations.(nameOutcome).(fSet).rs = rs;
    correlations.(nameOutcome).(fSet).p = p;
end
save('corr_RFprob_Time','correlations')

cd(startpath)
end