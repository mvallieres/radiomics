function displayBarResults_MetricBased_HN(pathExperiments,fSetNames,nameOutcome,metrics,titlePart,pathFig)
% - fSetNames: Cell of strings specifying the names of all feature sets
% - nameOutcome: String, specifying the name of a given specific outcome
% - titlePart: Put '' if no title addition.
% - pathFig: (optional). Full path to where figure is saved without
%   displaying it. Put '' for displaying the figure and not saving it to
%   'pathFig' (default).

startpath = pwd;

nExp = 1;

prefix = ['testResults',titlePart]; % So 5th argument is mandatory

nMetrics = numel(metrics);
if nargin < 6
    pathFig = '';
end

nFset = numel(fSetNames);
data = cell(1,nMetrics);
for m = 1:nMetrics
    data{m} = zeros(nExp,nFset);
end

% Obtaining data
cd(pathExperiments), pathExperiment = pwd;
for exp = 1:nExp
    for f = 1:nFset
        results = load([prefix,'_',fSetNames{f},'_',nameOutcome]); results = struct2cell(results); results = results{1};
        for m = 1:nMetrics
            data{m}(exp,f) = results.(metrics{m});
        end
    end
end

% Producing bar plot
scrsz = get(groot,'ScreenSize');
if isempty(pathFig)
    h = figure('Position',[1 7*scrsz(4)/10 19*scrsz(3)/40 7*scrsz(4)/10]);
else
    h = figure('visible','off','Position',[1 7*scrsz(4)/10 19*scrsz(3)/40 7*scrsz(4)/10]);
end
ind = strfind(nameOutcome,'Death');
if ~isempty(ind)
    nameOutcome(ind:ind+4) = [];
    nameOutcome = ['Survival',nameOutcome];
end

if ~isempty(strfind(fSetNames{1},'sign'))
    titleName = [nameOutcome,' -- ',titlePart,'sign'];
else
    titleName = [nameOutcome,' -- ',titlePart];
end

if exp > 1
    y = zeros(nMetrics,nFset);
    err = zeros(nMetrics,nFset);
    stdData = data;
    for m = 1:nMetrics
        data{m} = mean(data{m});
        stdData{m} = std(data{m});
        y(m,:) = data{m}(1,:);
        err(m,:) = stdData{m}(1,:);
    end
    model_series = y; model_error = err;
    bar(model_series,1);
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    set(gca,'XTicklabel','Modelo1|Modelo2|Modelo3')
    lh = legend('Serie1','Serie2','Serie3');
    set(lh,'Location','BestOutside','Orientation','horizontal')
    hold on;
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
          % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
          x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
          errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','LineWidth',3);
    end
else
    y = zeros(nMetrics,nFset);
    for m = 1:nMetrics
        y(m,:) = data{m}(1,:);
    end
    if nMetrics > 1
        bar(y,1)
    else
        bar([y;zeros(nFset-1,nFset)],1)
    end
end
set(gca,'XTickLabel',metrics);
set(gca,'FontSize',22)
axis([0.5 nMetrics+0.5 0 1])
legend(fSetNames)
title(titleName,'FontSize',26,'FontWeight','bold')
annotation(h,'line',[0.13 0.905],[0.515673981191221 0.515673981191221],...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'LineWidth',5,'LineStyle',':');

if ~isempty(pathFig)
    cd(pathFig)
    saveas(h,titleName,'fig')
end

cd(startpath)
end