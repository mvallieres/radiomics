function displayKaplanMeierRF_HNsupp(pathResults,nameOutcomes,varType,pathFig,thresholds)
% thresholds: (optional). Defines the thresholds of the risk stratification groups.
%             Defaut is one threshold of 0.5 for two risk groups.
% varType: Must either be 'radiomics_CT', 'radiomics_PET', 'radiomics_PETCT' or 'clinical'

if nargin < 5
    thresholds = 0.5;
else
    if numel(thresholds) > 2
        error('Max 3 curves (risk groups) are allowed for the moment')
    end
    thresholds = thresholds./100;
end
if numel(thresholds) == 1
    symbols = {'-b','--r'};
elseif numel(thresholds) == 2
    symbols = {'-g','--b','.-r'};
end

startpath = pwd;
nOutcomes = numel(nameOutcomes);

cd(pathResults), load('testing')
thresholds = [0,thresholds,1];
nCurves = numel(thresholds) - 1;
for o = 1:nOutcomes
    legendCell = cell(1,nCurves);
    nameOutcome = nameOutcomes{o};
    if strcmp(nameOutcome,'Locoregional')
        name = 'Recurrence-free probability';
    elseif strcmp(nameOutcome,'Distant')
        name = 'Metastasis-free probability';
    elseif strcmp(nameOutcome,'Death')
        name = 'Survival probability';
    elseif strcmp(nameOutcome,'DeathSign')
        name = 'Survival probability';
    end
    load(['testResultsRF',varType,'_',nameOutcome]) % results gets out of there
    prob = results.probResponse;
    if strcmp(nameOutcome,'DeathSign')
        outcome = testing.outcomes.Death;
        timeData =  testing.timeToEvents.Death;        
    else
        outcome = testing.outcomes.(nameOutcome);
        timeData =  testing.timeToEvents.(nameOutcome);
    end
    censData = 1 - outcome;
    if isempty(pathFig)
        hFig = figure;
    else
        hFig = figure('visible','off');
    end
    hold on, maxTT = 0; handles = []; tables1 = cell(1,nCurves); tables2 = cell(1,nCurves); X = cell(1,nCurves);
    for c = 1:nCurves
        time = timeData(prob >= thresholds(c) & prob < thresholds(c+1));
        cens = censData(prob >= thresholds(c) & prob < thresholds(c+1));
        [tableT,tableTT,t,T,xcg,ycg,~] = kmplot([time,cens],0.05,1,0); tables1{c} = tableT; tables2{c} = tableTT; X{c} = [time,cens];
        S = stairs(t,T,symbols{c},'LineWidth',2); handles = [handles,S];
        h = plot(xcg,ycg,'k+','MarkerSize',6,'LineWidth',1);
        maxT = max(t);
        if maxT > maxTT
            maxTT = maxT;
        end
        legendCell{c} = [num2str(thresholds(c)),' <= prob_{RF} < ',num2str(thresholds(c+1))];
    end
    if nCurves == 2
        pVal = getPval(tables1{1},tables1{2},tables2{1},tables2{2},X{1},X{2});
        annotation(hFig,'textbox',...
        [0.2635 0.300000001206285 0.22142856536167 0.0642857130794299],...
        'String',{['p - value = ',num2str(pVal)]},...
        'LineStyle','none','FitBoxToText','off');
    elseif nCurves == 3
        % Between low-risk and medium-risk curves (see bottom-left value)
        pVal = getPval(tables1{1},tables1{2},tables2{1},tables2{2},X{1},X{2});
        annotation(hFig,'textbox',...
        [0.2635 0.300000001206285 0.22142856536167 0.0642857130794299],...
        'String',{['p - value = ',num2str(pVal)]},...
        'LineStyle','none','FitBoxToText','off');
        
        % Between medium-risk and high-risk curves (see bottom-right value)
        pVal = getPval(tables1{2},tables1{3},tables2{2},tables2{3},X{2},X{3});
        annotation(hFig,'textbox',...
        [0.6135 0.300000001206285 0.22142856536167 0.0642857130794299],...
        'String',{['p - value = ',num2str(pVal)]},...
        'LineStyle','none','FitBoxToText','off');
    end
    axis([0 maxTT 0 1.05]), axis square
    ind = strfind(nameOutcome,'Death');
    if ~isempty(ind)
        nameOutcome(ind:ind+4) = [];
        nameOutcome = ['Survival',nameOutcome];
    end
    set(gca,'FontSize',12)
    titleName = ['Kaplan-Meier (SUPP): ',nameOutcome,'_{',varType,'}'];
    title(titleName,'FontSize',14)
    xlabel('Time (days)','FontSize',14)
    ylabel(name,'FontSize',14)
    if nCurves == 3
        legendCell = {'Low-risk group','Medium-risk group','High-risk group'};
    else
    end
    legend(handles,legendCell,'Location','SouthWest')
    grid ON
    hold off
    
    if ~isempty(pathFig)
        cd(pathFig)
        saveas(hFig,[titleName,'_Strat',num2str(nCurves)],'fig')
        cd(pathResults)
    end     
end

cd(startpath)
end