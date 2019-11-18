function estimateClinicalRF_Old(pathRF,nBoot,seed)

startpath = pwd;

% INITIALIZATIONS
cd(pathRF)
load('training'), load('testing')
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
fSetNames = fieldnames(training.textures.(nameOutcomes{1})); nFset = numel(fSetNames);
nClinical = size(training.clinical.table,2);
clinicNames = training.clinical.table.Properties.VariableNames';


% COMPUTATIONS
for o = 1:nOutcomes
    rng(seed), [trainSets,testSets] = buildBootSet(training.outcomes.(nameOutcomes{o}),nBoot);
    outcome = training.outcomes.(nameOutcomes{o});
    for f = 1:nFset
        text = training.textures.(nameOutcomes{o}).(fSetNames{f}); nText = size(text,2);
        maxAUC = 0;
        for c = 1:nClinical
            ind = combnk(1:nClinical,c); nComb = size(ind,1);
            for t = 1: nComb
                tableComb = [text,training.clinical.table(:,ind(t,:))];
                cat = logical([zeros(1,nText),training.clinical.categories(ind(t,:))]);
                meanAUC = 0;
                for n = 1:nBoot
                    Ttrain = tableComb(trainSets(:,n),:); Ttest = tableComb(testSets{n},:);
                    Ytrain = outcome(trainSets(:,n)); Ytest = outcome(testSets{n});
                    [RF] = applyEnsembleRF_table(Ttrain,Ytrain,cat);
                    [prob] = predictRF(Ttest,RF);
                    [aucBoot,~,~,~] = calcPerformMetrics(prob,Ytest,0.5);
                    meanAUC = meanAUC + aucBoot;
                end
                meanAUC = meanAUC/nBoot;
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