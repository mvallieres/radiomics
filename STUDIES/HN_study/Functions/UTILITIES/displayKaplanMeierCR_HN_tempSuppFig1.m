function displayKaplanMeierCR_HN_tempSuppFig1(pathResults,nameOutcomes,nameSets)
% thresholds: (optional). Defines the thresholds of the risk stratification groups.
%             Defaut is one threshold of 0.5 for two risk groups.

symbols = {'-b','--r'};

startpath = pwd;
nOutcomes = numel(nameOutcomes);

nCurves = 2;
cd(pathResults), load('testing')
for o = 1:nOutcomes
    nameOutcome = nameOutcomes{o}; fSet = nameSets{o};
    if strcmp(nameOutcome,'Locoregional')
        name = 'Recurrence-free probability';
    elseif strcmp(nameOutcome,'Distant')
        name = 'Metastasis-free probability';
    elseif strcmp(nameOutcome,'Death')
        name = 'Survival probability';
    elseif strcmp(nameOutcome,'DeathSign')
        name = 'Survival probability';
    end
    load(['testResultsCR_',fSet,'_',nameOutcome]) % results gets out of there
    resp = results.testData.response;
    outcome = results.testData.outcome;
    timeData =  results.testData.timeToEvent;
    censData = 1 - outcome;
    medianHR = results.model.medianHR;
    thresholds = [-Inf,medianHR,Inf];
    hFig =  subplot(1,nOutcomes,o); hold on, maxTT = 0; handles = []; tables1 = cell(1,nCurves); tables2 = cell(1,nCurves); X = cell(1,nCurves);
    for c = 1:nCurves
        time = timeData(resp >= thresholds(c) & resp < thresholds(c+1));
        cens = censData(resp >= thresholds(c) & resp < thresholds(c+1));
        [tableT,tableTT,t,T,xcg,ycg,~] = kmplot([time,cens],0.05,1,0); tables1{c} = tableT; tables2{c} = tableTT; X{c} = [time,cens];
        S = stairs(t,T,symbols{c},'LineWidth',2); handles = [handles,S];
        h = plot(xcg,ycg,'k+','MarkerSize',6,'LineWidth',1);
        maxT = max(t);
        if maxT > maxTT
            maxTT = maxT;
        end
    end

%     pVal = getPval(tables1{1},tables1{2},tables2{1},tables2{2},X{1},X{2});
%     annotation(hFig,'textbox',...
%     [0.2635 0.300000001206285 0.22142856536167 0.0642857130794299],...
%     'String',{['p - value = ',num2str(pVal)]},...
%     'LineStyle','none','FitBoxToText','off');
    
    legendCell = {'< medianHR','>= medianHR'};
    axis([0 maxTT 0 1.05]), axis square
    ind = strfind(nameOutcome,'Death');
    if ~isempty(ind)
        nameOutcome(ind:ind+4) = [];
        nameOutcome = ['Survival',nameOutcome];
    end
    set(gca,'FontSize',12)
    title(['Kaplan-Meier (CR): ',nameOutcome,'_{(',nameSets{o},')}'],'FontSize',14)
    xlabel('Time (days)','FontSize',14)
    ylabel(name,'FontSize',14)
    legend(handles,legendCell,'Location','SouthWest')
    grid ON
    hold off
end

cd(startpath)
end