function displayProbGraph_HN(pathRF,nameOutcomes,nameSets,thresholds,pathFig)
% nameOutcomes and nameSets are associated: one set per outcome
% - pathFig: (optional). Full path to where figure is saved without
%   displaying it. Put '' for displaying the figure and not saving it to
%   'pathFig' (default).

if nargin < 5
    pathFig = '';
end

startpath = pwd;
nOutcomes = numel(nameOutcomes);
symbols = {'ob','xr'}; % First entry for positive instances, second entry for negative instances

if isempty(pathFig)
    figure
else
    h = figure('visible','off');
end
cd(pathRF), load('testing')
nameThick = cell(1,nOutcomes*2+3); nameThick{1} = '';
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o}; fSet = nameSets{o};
    load(['testResultsRF_',fSet,'_',nameOutcome]) % results gets out of there
    probData = results.probResponse;
    if strcmp(nameOutcome,'DeathSign')
        outcome = testing.outcomes.Death;
    else
        outcome = testing.outcomes.(nameOutcome);
    end
    probPos = probData(outcome == 1); xPos = ones(numel(probPos),1)*o;
    probNeg = probData(outcome == 0); xNeg = ones(numel(probNeg),1)*o;
    plot(xPos,probPos,symbols{1},'LineWidth',6,'MarkerSize',18,'MarkerFaceColor',symbols{1}(end),'MarkerEdgeColor',symbols{1}(end))
    hold on
    plot(xNeg,probNeg,symbols{2},'LineWidth',6,'MarkerSize',20,'MarkerFaceColor',symbols{2}(end),'MarkerEdgeColor',symbols{2}(end))
    hold on
    nameThick{o*2} = ''; nameThick{o*2+1} = [nameOutcome,'_{',nameSets{o},'}'];
end
xThres = (0+0.05):0.05:(nOutcomes+1-0.05);
plot(xThres,ones(1,numel(xThres))*thresholds(1)/100,'--m','LineWidth',4), hold on
plot(xThres,ones(1,numel(xThres))*thresholds(2)/100,'--m','LineWidth',4)
nameThick{end} = '';
titleName = 'Probability of occurence of events: testing cohorts';
title(titleName,'FontSize',36)
ylabel('Random forest output probability (prob_{RF})','FontSize',24)
axis([0,nOutcomes+1,0,1])
legend('Status: Event occured','Status: Event did not occur')
set(gca,'xticklabel',nameThick)
set(gca,'FontSize',24)

if ~isempty(pathFig)
    cd(pathFig)
    saveas(h,titleName,'fig')
end 

cd(startpath)
end