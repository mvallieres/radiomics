function [prob] = predictRF_temp(X,Y,RF)

nInst = numel(Y);
trees = RF.Trees; numTrees = numel(trees);
predictions = zeros(nInst,1);
for t = 1:numTrees
    temp = trees{t}.predict(X); temp = str2double(temp);
    predictions = predictions + temp;
end
prob = predictions/numTrees;
end