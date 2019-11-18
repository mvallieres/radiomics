function estimateAllClinicalRF(pathRF,nSplit,nBoot,seed)

startpath = pwd;

% INITIALIZATIONS
cd(pathRF)
training = load('training'); training = struct2cell(training); training = training{1};
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
%fSetNames = fieldnames(training.textures.(nameOutcomes{1})); nFset = numel(fSetNames);
fSetNames = {'PET','CT','CT'}; nFset = 1;% Here, only 'PET' is tested for Loco, 'CT for 'Distant' and 'CT' for 'Death' 'f' was replacef by 'o' in three lines)
nClinical = size(training.clinical.table,2);
clinicNames = training.clinical.table.Properties.VariableNames';
testSplit = 1/3; % Proportion of test cases in the splits
comb =  {{[1,2,3],...      % Loco       : Age, Subtype, T_Stage
        [1,2,4],...        % Loco       : Age, Subtype, T_Stage, N_Stage
        [1,2,5],...        % Loco       : Age, Subtype, TNM_Stage
        [1,2,3,4]},...     % Loco       : Age, Subtype, T_Stage, N_Stage
        {[1,2,3],...       % Distant       : Age, Subtype, T_Stage
        [1,2,4],...        % Distant       : Age, Subtype, N_Stage
        [1,2,5],...        % Distant       : Age, Subtype, TNM_Stage
        [1,2,3,4]},...     % Distant       : Age, Subtype, T_Stage, N_Stage
        {[1,2,3],...       % Death       : Age, Subtype, T_Stage
        [1,2,4],...        % Death       : Age, Subtype, N_Stage
        [1,2,5],...        % Death       : Age, Subtype, TNM_Stage
        [1,2,3,4]}};       % Death       : Age, Subtype, T_Stage, N_Stage
    
nComb = [4,4,4];
testCost = 0.4:0.1:2.5; nCost = numel(testCost); % Emphasis factor on positive instances during random forest training
         
% COMPUTATIONS
for o = 1:nOutcomes
    outcome = training.outcomes.(nameOutcomes{o});
    rng(seed), [trainSets,testSets] = getSplits(outcome,nSplit,testSplit);
    for f = 1:nFset
        text = training.textures.(nameOutcomes{o}).(fSetNames{o}); nText = size(text,2);
        maxMetric = 0;
        for c = 1:nCost
            cost = testCost(c);
            for t = 1:nComb(o)
                tableComb = [text,training.clinical.table(:,comb{o}{t})];
                cat = logical([zeros(1,nText),training.clinical.categories(comb{o}{t})]);
                meanAUC = 0; meanSens = 0; meanSpec = 0;
                for s = 1:nSplit
                    Ttrain = tableComb(trainSets(:,s),:); Ttest = tableComb(testSets(:,s),:);
                    Ytrain = outcome(trainSets(:,s)); Ytest = outcome(testSets(:,s));
                    rng(seed), [RF] = trainRF_table(Ttrain,Ytrain,cat,nBoot,cost);
                    [prob] = predictRF(Ttest,RF);
                    [aucSplit,sensSplit,specSplit,~] = calcPerformMetrics(prob,Ytest,0.5);
                    meanAUC = meanAUC + aucSplit;
                    meanSens = meanSens + sensSplit;
                    meanSpec = meanSpec + specSplit;
                end
                meanAUC = meanAUC/nSplit; meanSens = meanSens/nSplit; meanSpec = meanSpec/nSplit;
                metric = 0.5*meanAUC + 0.5*(1-abs(meanSens - meanSpec));
                %metric = meanAUC;
                if metric >= maxMetric
                    var = clinicNames(comb{o}{t});
                    bestCost = cost;
                    indMax = comb{o}{t};
                    maxMetric = metric;
                end
            end
        end
        %fprintf('\n  * Optimal added clinical variables for "%s", "%s" are: ',nameOutcomes{o},fSetNames{f});
        %for v = 1:numel(var)-1
        %    fprintf([var{v},', '])
        %end
        %fprintf([var{end}])
        training.clinical.bestAdd.(nameOutcomes{o}).(fSetNames{o}) = indMax;
        training.cost.(nameOutcomes{o}).(fSetNames{o}) = bestCost;
    end
end
save('training','training')
fprintf('\n')

cd(startpath)
end