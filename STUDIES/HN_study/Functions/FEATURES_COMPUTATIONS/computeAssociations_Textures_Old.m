function computeAssociations_Textures_Old(pathUnivariate,pathFeatures,cohorts,outcomes,timeToEvent,featType,textType,textName,nPatient,scale_mat,algo_cell,Ng_mat)

% THIS FUNCTION WAS MODIFIED ON SEPTEMBER 29 TO COMPUTE FUSED TEXTURE
% RESULTS. PREVIOUS VERSION (WITHOUT FUSED) IS WITH COMMENTED LINES.


startpath = pwd;
cd(pathUnivariate)
fid = fopen('Texture_SignificanceProportion.txt','w');

nameOutcomes = fieldnames(outcomes.(cohorts{1})); nOutcomes = numel(nameOutcomes); nCohort = numel(cohorts);
nFeatType = numel(featType); nTextType = numel(textType);
nText = 0;
for t = 1:nTextType
    nText = nText + numel(textName{t});
end

% featureNames = {'Best PET texture - GTVp'; 'Best PET texture - GTVtot'; 'Best CT texture - GTVp'; 'Best CT texture - GTVtot'}; nFeature = numel(featureNames); 
featureNames = {'Best PET texture - GTVp'; 'Best PET texture - GTVtot'; 'Best CT texture - GTVp'; 'Best CT texture - GTVtot'; 'Best PET/CT texture - GTVp'; 'Best PET/CT texture - GTVtot'}; nFeature = numel(featureNames); 
rsMat_Binary = zeros(nFeature,nOutcomes); pMat_Binary = zeros(nFeature,nOutcomes); pCell_Binary = cell(nFeature,nOutcomes); stringCell_Binary = cell(nFeature,nOutcomes); nameTextCell_Binary = cell(nFeature,nOutcomes);
rsMat_Time = zeros(nFeature,nOutcomes); pMat_Time = zeros(nFeature,nOutcomes); pCell_Time = cell(nFeature,nOutcomes); stringCell_Time = cell(nFeature,nOutcomes); nameTextCell_Time = cell(nFeature,nOutcomes);
% scans = {'PT','CT'}; nScans = numel(scans); nParams = [numel(scale_mat)*numel(algo_cell)*numel(Ng_mat),numel(scale_mat)*numel(algo_cell)*numel(Ng_mat)]; paramSizes = {[numel(scale_mat),numel(algo_cell),numel(Ng_mat)],[numel(scale_mat),numel(algo_cell),numel(Ng_mat)]};
scans = {'PT','CT','PTCT'}; nScans = numel(scans); nParams = [numel(scale_mat)*numel(algo_cell)*numel(Ng_mat),numel(scale_mat)*numel(algo_cell)*numel(Ng_mat),numel(CTweight_mat)*numel(scale_mat)*numel(algo_cell)*numel(Ng_mat)]; paramSizes = {[numel(scale_mat),numel(algo_cell),numel(Ng_mat)],[numel(scale_mat),numel(algo_cell),numel(Ng_mat)],[numel(CTweight_mat),numel(scale_mat),numel(algo_cell),numel(Ng_mat)]};
nPat = 0;
for c = 1:nCohort
    nPat = nPat + nPatient.(cohorts{c});
end
i = 0;
for scan = 1:nScans
    nParam = nParams(scan);
    paramSize = paramSizes{scan};
    textParamName = cell(1,nText*nParam); count = 0;
    for t = 1:numel(textType)
        for tStr = 1:numel(textName{t})
            for p = 1:nParam
                count = count + 1;
                if numel(paramSize) == 3
                    [aa,bb,cc] = ind2sub(paramSize,p);
                    textParamName{count} = [textType{t},'/',textName{t}{tStr},' -- ','Scale=',num2str(scale_mat(aa)),',Quant.algo=',algo_cell{bb},',Ng=',num2str(Ng_mat(cc))];
                elseif numel(paramSize) == 4 % Remove this if fused is not used.
                    [aa,bb,cc,dd] = ind2sub(paramSize,p);
                    textParamName{count} = [textType{t},'/',textName{t}{tStr},' -- ','CTweight=',num2str(CTweight_mat(aa),'%.2f'),',Scale=',num2str(scale_mat(bb)),',Quant.algo=',algo_cell{cc},',Ng=',num2str(Ng_mat(dd))];
                end
            end
        end
    end
    for type = 1:nFeatType
        i = i + 1; textMat = zeros(nPat,nText*nParam); countPat = 0;
        outcomeMat = zeros(nPat,nOutcomes);
        timeMat = zeros(nPat,nOutcomes);
        cd(pathFeatures)
        for c = 1:nCohort
            cohort = cohorts{c};
            text = load(['text_',cohort,'_',scans{scan},'_',featType{type}]); text = struct2cell(text); text = text{1};
            patientVect = (countPat+1):(countPat+nPatient.(cohort)); countPat = countPat + nPatient.(cohort);
            count = 0;
            for t = 1:numel(textType)
                for tStr = 1:numel(textName{t})
                    for p = 1:nParam
                        count = count + 1;
                        if numel(paramSize) == 3
                            [aa,bb,cc] = ind2sub(paramSize,p);
                            textMat(patientVect,count) = text{aa,bb,cc}.(textType{t}).(textName{t}{tStr}).Data;
                        elseif numel(paramSize) == 4 % Remove this if fused is not used.
                            [aa,bb,cc,dd] = ind2sub(paramSize,p);
                            textMat(patientVect,count) = text{aa,bb,cc,dd}.(textType{t}).(textName{t}{tStr}).Data;
                        end
                    end
                end
            end
            for j = 1:nOutcomes
                outcomeMat(patientVect,j) = outcomes.(cohort).(nameOutcomes{j});
                timeMat(patientVect,j) = timeToEvent.(cohort).(nameOutcomes{j});
            end
        end
        cd(pathUnivariate)
        for j = 1:nOutcomes
            outcome = outcomeMat(:,j);
            time = timeMat(:,j);
            [rsTemp,pTemp] = corr(textMat,outcome,'type','Spearman','rows','pairwise'); temp = abs(rsTemp);
            [significance] = benjamini_hochberg(pTemp,0.10); proportion = sum(significance)/numel(significance);
            if proportion > 0
                starB = '*';
            else
                starB = '';
            end
            fprintf(fid,['Proportion of significant textures: ',scans{scan},', ',featType{type},', ',nameOutcomes{j},', Binary: ',num2str(proportion),'\n']);
            [~,indMax] = max(temp);
            rsMat_Binary(i,j) = rsTemp(indMax); pMat_Binary(i,j) = pTemp(indMax); nameTextCell_Binary{i,j} = textParamName{indMax};
            pCell_Binary{i,j} = pTemp; 
            [rsTemp,pTemp] = corr(textMat,time,'type','Spearman','rows','pairwise'); temp = abs(rsTemp);
            [significance] = benjamini_hochberg(pTemp,0.10); proportion = sum(significance)/numel(significance);
            if proportion > 0
                starT = '*';
            else
                starT = '';
            end
            fprintf(fid,['Proportion of significant textures: ',scans{scan},', ',featType{type},', ',nameOutcomes{j},', Time: ',num2str(proportion),'\n']);
            [~,indMax] = max(temp);
            rsMat_Time(i,j) = rsTemp(indMax); pMat_Time(i,j) = pTemp(indMax); nameTextCell_Time{i,j} = textParamName{indMax};
            pCell_Time{i,j} = pTemp; 
            if pMat_Binary(i,j) < 0.01
                stringCell_Binary{i,j} = [starB,'rs = ',num2str(rsMat_Binary(i,j),'%.2f'),', p = ',num2str(pMat_Binary(i,j),'%.2i'),starB];
            else
                stringCell_Binary{i,j} = [starB,'rs = ',num2str(rsMat_Binary(i,j),'%.2f'),', p = ',num2str(pMat_Binary(i,j),'%.2f'),starB];
            end
            if pMat_Time(i,j) < 0.01
                stringCell_Time{i,j} = [starT,'rs = ',num2str(rsMat_Time(i,j),'%.2f'),', p = ',num2str(pMat_Time(i,j),'%.2i'),starT];
            else
                stringCell_Time{i,j} = [starT,'rs = ',num2str(rsMat_Time(i,j),'%.2f'),', p = ',num2str(pMat_Time(i,j),'%.2f'),starT];
            end
        end
    end
end
cd(pathUnivariate), fclose(fid);
Locoregional = stringCell_Binary(:,1); Distant = stringCell_Binary(:,2); Death = stringCell_Binary(:,3);
results.tableCorr = table(Locoregional,Distant,Death,'RowNames',featureNames);
Locoregional = nameTextCell_Binary(:,1); Distant = nameTextCell_Binary(:,2); Death = nameTextCell_Binary(:,3);
results.tableTextName = table(Locoregional,Distant,Death,'RowNames',featureNames);
results.rsMat = rsMat_Binary; results.pMat = pMat_Binary; results.pCell = pCell_Binary;
save('texturesBest_UniV_Binary','results'), clear results
Locoregional = stringCell_Time(:,1); Distant = stringCell_Time(:,2); Death = stringCell_Time(:,3);
results.tableCorr = table(Locoregional,Distant,Death,'RowNames',featureNames);
Locoregional = nameTextCell_Time(:,1); Distant = nameTextCell_Time(:,2); Death = nameTextCell_Time(:,3);
results.tableTextName = table(Locoregional,Distant,Death,'RowNames',featureNames);
results.rsMat = rsMat_Time; results.pMat = pMat_Time; results.pCell = pCell_Time;
save('texturesBest_UniV_Time','results'), clear results

cd(startpath)
end