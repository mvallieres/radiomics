function estimateClinicalRF_Old2(pathRF,nSplit,nBoot,seed)

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
        maxAUC = 0;
        for c = 1:nClinical
            ind = combnk(1:nClinical,c); nComb = size(ind,1);
            for t = 1:nComb
                tableComb = [text,training.clinical.table(:,ind(t,:))];
                cat = logical([zeros(1,nText),training.clinical.categories(ind(t,:))]);
                meanAUC = 0;
                for s = 1:nSplit
                    Ttrain = tableComb(trainSets(:,s),:); Ttest = tableComb(testSets(:,s),:);
                    Ytrain = outcome(trainSets(:,s)); Ytest = outcome(testSets(:,s));
                    [RF] = trainRF_table(Ttrain,Ytrain,cat,nBoot);
                    [prob] = predictRF(Ttest,RF);
                    [aucSplit,~,~,~] = calcPerformMetrics(prob,Ytest,0.5);
                    meanAUC = meanAUC + aucSplit;
                end
                meanAUC = meanAUC/nSplit;
                if meanAUC >= maxAUC
                    var = clinicNames(ind(t,:));
                    indMax = ind(t,:);
                    maxAUC = meanAUC;
                end
            end
        end
        fprintf('\n  * Optimal added clinical variables for "%s", "%s" are: ',nameOutcomes{o},fSetNames{f});
        for v = 1:numel(var)-1
            fprintf([var{v},', '])
        end
        fprintf([var{end}])
        testing.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{f}) = indMax;
    end
end
save('testing','testing')
fprintf('\n')

cd(startpath)
end