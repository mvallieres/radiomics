function testFinalModels_HN_OldExp(pathExperiments,nExp,fSetNames,paramCells,nameCells)
% - paramCells: One cell of parameters for each fSet
% - nameCells: One cell of names for each fSet

startpath = pwd;
nFset = numel(fSetNames);
for exp = 1:nExp
    cd(fullfile(pathExperiments,['Experiment',num2str(exp)])), pathExperiment = pwd;
    load('training'), load('testing')
    nameOutcomes = fieldnames(training); nOutcomes = numel(nameOutcomes); 
    cd('FINAL_MODELS'), pathFinalModels = pwd;
    for o = 1:nOutcomes
        for f = 1:nFset
            results = []; results = struct;
            cd(fullfile(pathFinalModels,nameOutcomes{o},fSetNames{f}))
            load('finalModel'), load('coeff'), load('response'), model = struct;
            
            % OBTAINING MODEL CHARACTERISTIC
            order = numel(finalModel.Name);
            paramCell = paramCells{f}; nameCell = nameCells{f};
            nParam = numel(paramCell); scans = cell(1,order); 
            for n = 1:order
                name = finalModel.Name{n};
                indA = strfind(name,'--'); indB = strfind(name,'-');
                if ~isempty(indA) % This is a texture
                    model.(['Feature',num2str(n)]).isTexture = true;
                    model.(['Feature',num2str(n)]).textType = name((indA(1)+2):(indB(3)-1));
                    model.(['Feature',num2str(n)]).textName = name((indB(3)+1):end);
                    indScan = strfind(name,'(');
                    scans{n} = name(1:indScan-1);
                    model.(['Feature',num2str(n)]).scan = scans{n};
                    indMat = zeros(1,nParam); sizeParam = zeros(1,nParam);
                    for p = 1:nParam
                        paramAll = paramCell{p}; nSpecParam = numel(paramAll); sizeParam(p) = nSpecParam;
                        nameParam = nameCell{p};
                        indP = strfind(name,nameParam);
                        for start = indP:numel(name)
                            if strcmp(name(start),'=')
                                break
                            end
                        end
                        for final = indP:numel(name)
                            if strcmp(name(final),',') | strcmp(name(final),')')
                                break
                            end
                        end
                        if iscell(paramCell{p})
                            param = name((start+1):(final-1));
                            for nsp = 1:nSpecParam
                                if strcmp(param,paramAll{nsp})
                                    indMat(p) = nsp;
                                    break
                                end
                            end
                        else
                            param = str2num(name((start+1):(final-1)));
                            for nsp = 1:nSpecParam
                                if param > (paramAll(nsp)-0.05) & param < (paramAll(nsp)+0.05)
                                    indMat(p) = nsp;
                                    break
                                end
                            end
                        end
                    end
                    if numel(sizeParam) == 2
                        model.(['Feature',num2str(n)]).indText = sub2ind(sizeParam,indMat(1),indMat(2));
                    elseif numel(sizeParam) == 3
                        model.(['Feature',num2str(n)]).indText = sub2ind(sizeParam,indMat(1),indMat(2),indMat(3));
                    elseif numel(sizeParam) == 4
                        model.(['Feature',num2str(n)]).indText = sub2ind(sizeParam,indMat(1),indMat(2),indMat(3),indMat(4));
                    elseif numel(sizeParam) == 5
                        model.(['Feature',num2str(n)]).indText = sub2ind(sizeParam,indMat(1),indMat(2),indMat(3),indMat(4),indMat(5));
                    elseif numel(sizeParam) == 6
                        model.(['Feature',num2str(n)]).indText = sub2ind(sizeParam,indMat(1),indMat(2),indMat(3),indMat(4),indMat(5),indMat(6));
                    end
                    model.(['Feature',num2str(n)]).FullName = name;
                else % This is a non-texture
                    scans{n} = 'nonText';
                    model.(['Feature',num2str(n)]).isTexture = false;
                    model.(['Feature',num2str(n)]).nonText = name;
                end
            end
            model.coeff = coeff;
            model.order = order;
            results.model = model;
            results.trainData.data = finalModel.Data; results.trainData.response = response; results.trainData.outcome = training.(nameOutcomes{o}).outcome;
            
            % TESTING MODEL
            text = testing.(nameOutcomes{o}).text; nonText = testing.(nameOutcomes{o}).nonText;
            outcome = testing.(nameOutcomes{o}).outcome;
            resp = zeros(numel(outcome),1);
            data = zeros(numel(outcome),order);
            for n = 1:order
                coeff = model.coeff(n);
                feat = ['Feature',num2str(n)];
                scan = scans{n};
                if model.(feat).isTexture
                    data(:,n) = text.(fSetNames{f}).(scan){model.(feat).indText}.(model.(feat).textType).(model.(feat).textName).Data;
                else
                    data(:,n) = nonText.(model.(feat).nonText).Data;
                end
                resp = resp + coeff*data(:,n);
            end
            resp = resp + model.coeff(end);
            results.testData.data = data; results.testData.response = resp; results.testData.outcome = outcome;
            [AUC,sensitivity,specificity,accuracy] = calcPerformMetrics(resp,outcome,0);
            results.AUC = AUC;
            results.Sensitivity = sensitivity;
            results.Specificity = specificity;
            results.Accuracy = accuracy;
            cd(pathExperiment), save(['testResultsModels_',fSetNames{f},'_',nameOutcomes{o}],'results')
        end
    end
end

cd(startpath)
end