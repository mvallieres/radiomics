function [RF] = applyEnsembleRF_table(table,Y,categories,cost)


% INITIALIZATION
nNeg = sum(~Y); % Number of negative instances
nPos = sum(Y); % Number of positive instances
indPos = find(Y); % Indexes of positive instances
indNeg = find(~Y); % Indexes of negative instances
nSub = round(nNeg/nPos); % Number of subsets


% RANDOM NUMBER GENERATOR SEED
if ~RandStream.getGlobalStream.Seed
    rng('shuffle')
end

if nargin < 4
    cost = 1;
end

% FIND THE NUMBER OF NEGATIVE INSTANCES IN EACH SUBSET
nNegS = nNeg / nSub; nSup = ceil(nNegS); nInf = floor(nNegS);
if nSup ~= nInf
    nSubInf = nSub - 1; nSubSup = 1; total = nSubInf*nInf + nSubSup*nSup;
    while total ~= nNeg
        nSubInf = nSubInf - 1; nSubSup = nSubSup + 1;
        total = nSubInf*nInf + nSubSup*nSup;
    end
    nNegS = [repmat(nInf,[1,nSubInf]),repmat(nSup,[1,nSubSup])];
else % The number of negative instances in all partitions will be the same
    nNegS = repmat(nSup,[1,nSub]);
end


% FINDING THE INDEXES OF NEGATIVE INSTANCES IN EACH SUBSET
indNegSub = cell(1,nSub);
for i = 1:nSub-1
    indNegSub{i} = zeros(nNegS(i),1);
    indTemp = ceil(numel(indNeg)*rand(nNegS(i),1));
    indTemp = unique(indTemp);
    total = numel(indTemp);
    while total ~= nNegS(i)
        indMore = ceil(numel(indNeg)*rand(nNegS(i)-total,1));
        indTemp = [indTemp;unique(indMore)];
        indTemp = unique(indTemp);
        total = numel(indTemp);
    end
    indNegSub{i}(:) = indNeg(indTemp);
    indNeg(indTemp) = [];
end
indNegSub{end} = indNeg;


% COMPUTING A RANDOM FOREST FOR EACH SUBSET AND APPENDING
ind = [indNegSub{1};indPos];
Xtemp = table(ind,:);
Ytemp = Y(ind);
[Xtemp,Ytemp] = shufflePartition(Xtemp,Ytemp);
RF = TreeBagger(1,Xtemp,Ytemp,'CategoricalPredictors',categories,'SampleWithReplacement','on','Cost',[0,1/cost;1,0]);
for i = 2:nSub
    ind = [indNegSub{i};indPos];
    Xtemp = table(ind,:);
    Ytemp = Y(ind);
    [Xtemp,Ytemp] = shufflePartition(Xtemp,Ytemp);
    RFtemp = TreeBagger(1,Xtemp,Ytemp,'CategoricalPredictors',categories,'SampleWithReplacement','on','Cost',[0,1/cost;1,0]);
    RF = append(RF,RFtemp);
end

end