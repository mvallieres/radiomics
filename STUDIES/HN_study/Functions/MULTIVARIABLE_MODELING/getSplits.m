function [trainSets,testSets] = getSplits(outcome,nSplit,testSplit)

indNeg = find(outcome == 0); nNeg = numel(indNeg);
indPos = find(outcome == 1); nPos = numel(indPos);
nPosTest = round(testSplit*nPos); % Number of positive instances in the test sets
nNegTest = round(testSplit*nNeg); % Number of negative instances in the test sets
nInst = numel(outcome);
nTest = nPosTest + nNegTest;
nTrain = nInst - nTest;


trainSets = zeros(nTrain,nSplit);
testSets = zeros(nTest,nSplit);
for s = 1:nSplit
    indPosTest = indPos(datasample(1:nPos,nPosTest,'Replace',false)');
    indNegTest = indNeg(datasample(1:nNeg,nNegTest,'Replace',false)');
    indTest = [indPosTest;indNegTest]; indPerm = randperm(nTest)'; indTest = indTest(indPerm);
    indTrain = (1:nInst)'; indTrain(indTest) = []; indPerm = randperm(nTrain)'; indTrain = indTrain(indPerm);
    trainSets(:,s) = indTrain; testSets(:,s) = indTest;
end

end