function [indClinic,bestCost] = estimateClinicalRF(pathRF,nameOutcome,testComb,testCost,nSplit,testSplit,nBoot,seed,output)

% output: optional. 

if nargin < 9
    output = 1; % By default
end

startpath = pwd;
tStart = tic;

% INITIALIZATIONS
cd(pathRF)
rng(seed), seeds = ceil(1000000000*rand(2,1));
training = load('training'); training = struct2cell(training); training = training{1};
nComb = numel(testComb);
nCost = numel(testCost);
nTest = nComb *nCost;       

% COMPUTATIONS
if strcmp(nameOutcome,'DeathSign')
    outcome = training.outcomes.Death;
else
    outcome = training.outcomes.(nameOutcome);
end
rng(seeds(1)), [trainSets,testSets] = getSplits(outcome,nSplit,testSplit);
maxMetric = 0; count = 0;
for c = 1:nCost
    cost = testCost(c);
    for t = 1:nComb
        tableComb = training.clinical.table(:,testComb{t});
        cat = logical(training.clinical.categories(testComb{t}));
        meanAUC = 0; meanSens = 0; meanSpec = 0; count = count + 1;
        if output
            tic, fprintf('\n--> ESTIMATING RF PERFORMANCE USING %u SPLITS FOR TEST %u OF %u: %s ... ',nSplit,count,nTest,nameOutcome)
        end
        for s = 1:nSplit
            Ttrain = tableComb(trainSets(:,s),:); Ttest = tableComb(testSets(:,s),:);
            Ytrain = outcome(trainSets(:,s)); Ytest = outcome(testSets(:,s));
            rng(seeds(2)), [RF] = trainRF_table(Ttrain,Ytrain,cat,nBoot,cost);
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
            bestCost = cost;
            indClinic = testComb{t};
            maxMetric = metric;
        end
        if output
            fprintf('DONE!\n'), toc
        end
    end
end
if output
    time = toc(tStart);
    fprintf('\n\n\n-------------------------------------\n')
    fprintf('TOTAL TIME: %f seconds',time)
    fprintf('\n-------------------------------------\n')
end

cd(startpath)
end