function displayProbGraph(pathRF,nameOutcomes,nameSets,thesholds)
% nameOutcomes and nameSets are associated: one set per outcome

startpath = pwd;
nOutcomes = numel(nameOutcomes);
symbols = {'ob','xr'}; % First entry for positive instances, second entry for negative instances

figure
cd(pathRF), load('testing')
nameThick = cell(1,nOutcomes*2+3); nameThick{1} = '';
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o}; fSet = nameSets{o};
    %load(['testResultsRF_',fSet,'_',nameOutcome]) % results gets out of there
    load(['testResults_',fSet,'_',nameOutcome]) % results gets out of there
    probData = results.probResponse;
    outcome = testing.outcomes.(nameOutcome);
    probPos = probData(outcome == 1); xPos = ones(numel(probPos),1)*o;
    probNeg = probData(outcome == 0); xNeg = ones(numel(probNeg),1)*o;
    plot(xPos,probPos,symbols{1},'LineWidth',6,'MarkerSize',18,'MarkerFaceColor',symbols{1}(end),'MarkerEdgeColor',symbols{1}(end))
    hold on
    plot(xNeg,probNeg,symbols{2},'LineWidth',6,'MarkerSize',20,'MarkerFaceColor',symbols{2}(end),'MarkerEdgeColor',symbols{2}(end))
    hold on
    nameThick{o*2} = ''; nameThick{o*2+1} = [nameOutcome,'_{',nameSets{o},'+clinic}'];
end
xThres = (0+0.05):0.05:(nOutcomes+1-0.05);
plot(xThres,ones(1,numel(xThres))*thesholds(1)/100,'--m','LineWidth',4), hold on
plot(xThres,ones(1,numel(xThres))*thesholds(2)/100,'--m','LineWidth',4)
nameThick{end} = '';
title('Probability of occurence of events: testing cohorts','FontSize',36)
ylabel(['Random forest output probability'],'FontSize',24)
axis([0,nOutcomes+1,0,1])
legend('Status: Event occured','Status: Event did not occur')
set(gca,'xticklabel',nameThick)
set(gca,'FontSize',24)

cd(startpath)
end