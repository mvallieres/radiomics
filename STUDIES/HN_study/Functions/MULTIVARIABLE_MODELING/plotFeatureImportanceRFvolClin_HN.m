function plotFeatureImportanceRFvolClin_HN(pathRF)

startpath = pwd;

cd(pathRF), load('training'), load('testingVariableImportance')
nameOutcomes = fieldnames(training.outcomes); nOutcomes = numel(nameOutcomes);
% fSetNames = fieldnames(training.textures.(nameOutcomes{1})); nFset = numel(fSetNames);
fSetNames = {'volClinical'}; nFset = numel(fSetNames);

for o = 1:nOutcomes
    for f = 1:nFset
        percentAUCdecrease = variableImportance.(nameOutcomes{o}).(fSetNames{f}).percentAUCdecrease;
        varNames = variableImportance.(nameOutcomes{o}).(fSetNames{f}).varNames; nVar = numel(varNames);
        figure
        barh(1:nVar,percentAUCdecrease*100)
        for i = 1:nVar
            ind = strfind(varNames{i},'_'); nInd = numel(ind);
            if ~isempty(ind)
                for n = 1:nInd
                    varNames{i} = [varNames{i}(1:(ind(n)-1+(n-1))),'\',varNames{i}((ind(n)+(n-1)):end)];
                end
            end
        end
        set(gca,'yticklabel',varNames)
        nameOutcome = nameOutcomes{o};
        ind = strfind(nameOutcome,'Death');
        if ~isempty(ind)
            nameOutcome(ind:ind+4) = [];
            nameOutcome = ['Survival',nameOutcome];
        end
        title(['RANDOM PERMUTATIONS (RF):',nameOutcome,' -- ',fSetNames{f}])
        xlabel('Percent AUC decrease')
    end
end

cd(startpath)
end