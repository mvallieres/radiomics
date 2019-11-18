function [weights] = findFusionWeights(pathWORK,nPatient,outcomes,roiNumb,wavelet,nBoot)
% TESTING ISSAM'S METHOD FOR PET/CT fusion
%
% - pathWORK: HN WORKSPACE where the 'DATA' folder is located
% - nPatient: Number of patients
% - outcomes: Structure defining all outcomes.
% - roiNumb: 1st column for PET, 2nd column for CT
% - wavelet: String for the name of the wavelet ('sym8')
% - nBoot: (optional). 
% - weights: PET-HLL, CT-HLL, PET-LHL, CT-LHL, PET-HHL, CT-HHL, PET-LLH,
%            CT-LLH, PET-HLH, CT-HLH, PET-LHH, CT-LHH, PET-HHH, CT-HHH
%            First column is first outcome in outcomes, second column is second
%            outcome, etc.

startpath = pwd;
cd([pathWORK,'/DATA'])

% INTIALIZATION
bands = 2:8;
thresh = 0.01;
nBands = length(bands);
nameOutcomes = fieldnames(outcomes);
nOutcomes = length(nameOutcomes);
weights = zeros(nBands*2,nOutcomes);
v0 = cell(1,nOutcomes); % Cell of matrices of negative instances, for each outcome. Each column corresponds to a sub-band of a scan (total = 2 scans * 7 sub-bands)
v1 = cell(1,nOutcomes); % Cell of matrices of positive instances, for each outcome. Each column corresponds to a sub-band of a scan (total = 2 scans * 7 sub-bands)

% LOADING AND PRE-PROCESSING ALL VOLUMES
wdecPET = cell(1,nPatient);
wdecCT = cell(1,nPatient);
for i = 1:nPatient
    load(['Patient',num2str(i),'_PET'])
    box = getROIbox(sData,roiNumb(i,1));
    szPET = size(box);
    mask = sData{2}.scan.contour(roiNumb(i,1)).boxMask;
    box = box - min(box(mask==1));
    box = box/max(box(mask==1))*255 + 1; 
    box = box.*mask; % Setting values outside the ROI to 0, anything else between 1 and 256
    wdecPET{i} = wavedec3(box,1,wavelet);
    
    load(['Patient',num2str(i),'_CT'])
    box = getROIbox(sData,roiNumb(i,2));
    mask = sData{2}.scan.contour(roiNumb(i,2)).boxMask;
    tempBox = zeros(szPET);
    tempMask = zeros(szPET);
    for s = 1:szPET(3)
        tempBox(:,:,s) = imresize(box(:,:,s),[szPET(1),szPET(2)],'Method','cubic','Antialiasing',true);
        tempMask(:,:,s) = imresize(mask(:,:,s),[szPET(1),szPET(2)],'nearest');
    end
    box = tempBox; mask = tempMask;
    box = box - min(box(mask==1));
    box = box/max(box(mask==1))*255 + 1; 
    box = box.*mask; % Setting values outside the ROI to 0, anything else between 1 and 256
    wdecCT{i} = wavedec3(box,1,wavelet);
end

% COMPUTING WEIGHTS FOR ALL OUTCOMES
for o = 1:nOutcomes
    outcome = outcomes.(nameOutcomes{o});
    if nargin == 6
        [bootSam,~] = buildBootSet(outcome,nBoot,'adjust'); % 'trainSets' is a matrix, 'testSets' is a cell 
    else
        nBoot = 1;
        bootSam = (1:nPatient)';
    end
    for n = 1:nBoot
        v0{o} = []; v1{o} = [];
        for i = 1:nPatient
            tempData = zeros(numel(wdecPET{bootSam(i,n)}.dec{1}),nBands*2);
            for b = 1:nBands
                tempData(:,(b-1)*2 + 1) = wdecPET{bootSam(i,n)}.dec{bands(b)}(:);
                tempData(:,b*2) = wdecCT{bootSam(i,n)}.dec{bands(b)}(:);
            end
            ind = [];
            for b = 1:nBands*2
                ind = [ind;find(abs(tempData(:,b))<thresh)];
            end
            ind = unique(ind);
            tempData(ind,:) = [];
            if outcome(bootSam(i,n))
                v1{o} = [v1{o};tempData];
            else
                v0{o} = [v0{o};tempData];
            end
        end
        v1M = mean(v1{o});
        v0M = mean(v0{o});
        for j = 1:nBands*2
            v1{o}(:,j) = v1{o}(:,j) - v1M(j);
            v0{o}(:,j) = v0{o}(:,j) - v1M(j);
        end
        S = 0.5*(v0{o}'*v0{o}) + 0.5*(v1{o}'*v1{o});
        weights(:,o) = weights(:,o) + ((v1M - v0M)*S^(-1))';
    end
end
weights = weights/nBoot;

% RENORMALIZING WEIGHTS (we want the sum of the weights to equal 7)
for o = 1:nOutcomes
    weights(:,o) = weights(:,o) - min(weights(:,o));
    weights(:,o) = weights(:,o) + 0.01*max(weights(:,o));
    weights(:,o) = weights(:,o)*(nBands/(sum(weights(:,o))));
end

cd(startpath)
end