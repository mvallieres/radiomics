function computeAssociations_Textures(pathUnivariate,pathFeatures,cohorts,outcomes,featType,textType,textName,nPatient,scale_mat,algo_cell,Ng_mat)

startpath = pwd;
cd(pathUnivariate)
fid = fopen('Texture_SignificanceProportion.txt','w');

nameOutcomes = fieldnames(outcomes.(cohorts{1})); nOutcomes = numel(nameOutcomes); nCohort = numel(cohorts);
nFeatType = numel(featType); nTextType = numel(textType);
nText = 0;
for t = 1:nTextType
    nText = nText + numel(textName{t});
end

featureNames = {'Best PET texture - GTVp'; 'Best PET texture - GTVtot'; 'Best CT texture - GTVp'; 'Best CT texture - GTVtot'}; nFeature = numel(featureNames); 
rsMat_Binary = zeros(nFeature,nOutcomes); pMat_Binary = zeros(nFeature,nOutcomes); pCell_Binary = cell(nFeature,nOutcomes); stringCell_Binary = cell(nFeature,nOutcomes); nameTextCell_Binary = cell(nFeature,nOutcomes);
scans = {'PT','CT'}; nScans = numel(scans); nParams = [numel(scale_mat)*numel(algo_cell)*numel(Ng_mat),numel(scale_mat)*numel(algo_cell)*numel(Ng_mat)]; paramSizes = {[numel(scale_mat),numel(algo_cell),numel(Ng_mat)],[numel(scale_mat),numel(algo_cell),numel(Ng_mat)]};
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
                [aa,bb,cc] = ind2sub(paramSize,p);
                textParamName{count} = [textType{t},'/',textName{t}{tStr},' -- ','Scale=',num2str(scale_mat(aa)),',Quant.algo=',algo_cell{bb},',Ng=',num2str(Ng_mat(cc))];
            end
        end
    end
    for type = 1:nFeatType
        i = i + 1; textMat = zeros(nPat,nText*nParam); countPat = 0;
        outcomeMat = zeros(nPat,nOutcomes);
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
                        [aa,bb,cc] = ind2sub(paramSize,p);
                        textMat(patientVect,count) = text{aa,bb,cc}.(textType{t}).(textName{t}{tStr}).Data;
                    end
                end
            end
            for j = 1:nOutcomes
                outcomeMat(patientVect,j) = outcomes.(cohort).(nameOutcomes{j});
            end
        end
        cd(pathUnivariate)
        for j = 1:nOutcomes
            outcome = outcomeMat(:,j);
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
            if pMat_Binary(i,j) < 0.01
                stringCell_Binary{i,j} = [starB,'rs = ',num2str(rsMat_Binary(i,j),'%.2f'),', p = ',num2str(pMat_Binary(i,j),'%.2i'),starB];
            else
                stringCell_Binary{i,j} = [starB,'rs = ',num2str(rsMat_Binary(i,j),'%.2f'),', p = ',num2str(pMat_Binary(i,j),'%.2f'),starB];
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

cd(startpath)
end