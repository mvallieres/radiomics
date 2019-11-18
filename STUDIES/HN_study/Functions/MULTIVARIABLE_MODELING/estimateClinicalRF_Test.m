function estimateClinicalRF_Test(pathRF,nSplit,nBoot,seed)

startpath = pwd;

% INITIALIZATIONS
cd(pathRF)
training = load('training'); training = struct2cell(training); training = training{1};
testing = load('testing'); testing = struct2cell(testing); testing = testing{1};
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
fSetNames = fieldnames(training.textures.(nameOutcomes{1})); nFset = numel(fSetNames);
nClinical = size(training.clinical.table,2);
clinicNames = training.clinical.table.Properties.VariableNames';
testSplit = 1/3; % Proportion of test cases in the splits

         
% COMPUTATIONS
for o = 1:nOutcomes
    outcome = training.outcomes.(nameOutcomes{o});
    rng(seed), [trainSets,testSets] = getSplits(outcome,nSplit,testSplit);
    for f = 1:nFset
        text = training.textures.(nameOutcomes{o}).(fSetNames{f}); nText = size(text,2);
        indLeft = 1:nClinical;
        indChosen = [];
        for t = 1:nClinical
            maxMetric = 0;
            for i = 1:(nClinical-t+1)
                tableComb = [text,training.clinical.table(:,[indChosen,indLeft(i)])];
                cat = logical([zeros(1,nText),training.clinical.categories([indChosen,indLeft(i)])]);
                meanAUC = 0; meanSens = 0; meanSpec = 0;
                for s = 1:nSplit
                    Ttrain = tableComb(trainSets(:,s),:); Ttest = tableComb(testSets(:,s),:);
                    Ytrain = outcome(trainSets(:,s)); Ytest = outcome(testSets(:,s));
                    rng(seed), [RF] = trainRF_table(Ttrain,Ytrain,cat,nBoot);
                    [prob] = predictRF(Ttest,RF);
                    [aucSplit,sensSplit,specSplit,~] = calcPerformMetrics(prob,Ytest,0.5);
                    meanAUC = meanAUC + aucSplit;
                    meanSens = meanSens + sensSplit;
                    meanSpec = meanSpec + specSplit;
                end
                meanAUC = meanAUC/nSplit; meanSens = meanSens/nSplit; meanSpec = meanSpec/nSplit;
                %metric = 0.5*meanAUC + 0.5*(1-abs(meanSens - meanSpec));
                metric = meanAUC;
                if metric >= maxMetric
                    choice = indLeft(i);
                    indBye = i;
                    maxMetric = metric;
                end
            end
            indChosen = [indChosen,choice];
            var = clinicNames(indChosen);
            fprintf('\n  * Optimal added clinical variables for order %u, "%s", "%s" are: ',t,nameOutcomes{o},fSetNames{f});
            for v = 1:numel(var)-1
                fprintf([var{v},', '])
            end
            fprintf([var{end}])
            testing.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{f}).(['Order',num2str(t)]).Indexes = indChosen;
            testing.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{f}).(['Order',num2str(t)]).AUC = maxMetric;
            indLeft(indBye) = [];
        end
    end
end
save('testing','testing')
fprintf('\n')

cd(startpath)
end