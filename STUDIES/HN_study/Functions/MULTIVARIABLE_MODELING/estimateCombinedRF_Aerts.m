function estimateCombinedRF_Aerts(pathRF,nSplit,nBoot,seed)

startpath = pwd;

% INITIALIZATIONS
cd(pathRF), load('training')
testSplit = 1/3; % Proportion of test cases in the splits
testCost = 0.5:0.1:2; % Emphasis factor on positive instances during random forest training
comb = {[1,2,3],[1,2,4],[1,2,5],[1,2,3,4]};

% COMPUTATIONS
[indClinic,bestCost] = estimateCombinedRF(pathRF,'DeathSign','CT',comb,testCost,nSplit,testSplit,nBoot,seed,0);
training.clinical.bestAdd.DeathSign.CT = indClinic;
training.cost.DeathSign.CT = bestCost;

cd(pathRF), save('training','training')
cd(startpath)
end