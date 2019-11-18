function plotFeatureImportanceRF_HN(pathRF,pathFig)
% - pathFig: (optional). Full path to where figure is saved without
%   displaying it. Put '' for displaying the figure and not saving it to
%   'pathFig' (default).


if nargin < 2
    pathFig = '';
end

% VARYING PARAMETERS, DEPENDENT ON FINAL CHOICE OF BEST RESULTS (one Rad/Vol set for each outcome + clinical variables)
fSetNames = {'PETCTclinic','CTclinic','clinic'}; % For Locoregional, Distant and Death outcome, respectively
nameOutcomes = {'Locoregional','Distant','Death'};
nOutcomes = numel(nameOutcomes);


startpath = pwd;
cd(pathRF), load('training'), load('testingVariableImportance')

for o = 1:nOutcomes
    percentAUCchange = variableImportance.(nameOutcomes{o}).(fSetNames{o}).percentAUCchange;
    varNames = variableImportance.(nameOutcomes{o}).(fSetNames{o}).varNames; nVar = numel(varNames);
    if isempty(pathFig)
        figure
    else
        h = figure('visible','off');
    end
    barh(1:nVar,percentAUCchange*100)
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
    titleName = ['RANDOM PERMUTATIONS (RF):',nameOutcome,' -- ',fSetNames{o}];
    title(titleName)
    xlabel('Percent AUC change')
    
    if ~isempty(pathFig)
        cd(pathFig)
        saveas(h,titleName,'fig')
        cd(pathRF)
    end 
end

cd(startpath)
end