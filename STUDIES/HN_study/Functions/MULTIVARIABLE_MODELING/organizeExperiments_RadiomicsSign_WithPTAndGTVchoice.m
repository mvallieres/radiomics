function organizeExperiments_RadiomicsSign_WithPTAndGTVchoice(pathFeatures,pathLR,pathCR,pathRF,outcomes)

startpath = pwd;

% ALL FOR GTVtot

cohorts = fieldnames(outcomes);
nCohort = numel(cohorts);
% for c = 1:nCohort
%     outcomes.(cohorts{c}) = rmfield(outcomes.(cohorts{c}),{'Locoregional','Distant'});
% end
nameOutcomes = fieldnames(outcomes.(cohorts{1})); nOutcomes = numel(nameOutcomes);
fSetNames = {'PET','CT'}; % Only those are tested for the radiomic signature.
nFset = numel(fSetNames);


% FOR LR FOLDER
cd(pathLR), load('training'), load('testing'), load('permTrain'), load('permTest')
cd(pathFeatures), signTrain = struct; signTest = struct;
train1 = load('AertsSign_HGJ_PT_GTVtot'); train1 = struct2cell(train1); train1 = train1{1};
train2 = load('AertsSign_CHUS_PT_GTVtot'); train2 = struct2cell(train2); train2 = train2{1};
test1 = load('AertsSign_HMR_PT_GTVtot'); test1 = struct2cell(test1); test1 = test1{1};
test2 = load('AertsSign_CHUM_PT_GTVtot'); test2 = struct2cell(test2); test2 = test2{1};
signTrain.PET = [train1.Data;train2.Data]; signTrain.PET = signTrain.PET(permTrain,:);
signTest.PET = [test1.Data;test2.Data]; signTest.PET = signTest.PET(permTest,:);
train1 = load('AertsSign_HGJ_CT_GTVtot'); train1 = struct2cell(train1); train1 = train1{1};
train2 = load('AertsSign_CHUS_CT_GTVtot'); train2 = struct2cell(train2); train2 = train2{1};
test1 = load('AertsSign_HMR_CT_GTVtot'); test1 = struct2cell(test1); test1 = test1{1};
test2 = load('AertsSign_CHUM_CT_GTVtot'); test2 = struct2cell(test2); test2 = test2{1};
signTrain.CT = [train1.Data;train2.Data]; signTrain.CT = signTrain.CT(permTrain,:);
signTest.CT = [test1.Data;test2.Data]; signTest.CT = signTest.CT(permTest,:);
for o = 1:nOutcomes
    training.(nameOutcomes{o}).sign = signTrain;
    testing.(nameOutcomes{o}).sign = signTest;
end
cd(pathLR), save('training','training'), save('testing','testing')
clear training testing


% FOR CR FOLDER
cd(pathCR), load('training'), load('testing'), load('permTrain'), load('permTest')
cd(pathFeatures), signTrain = struct; signTest = struct;
train1 = load('AertsSign_HGJ_PT_GTVtot'); train1 = struct2cell(train1); train1 = train1{1};
train2 = load('AertsSign_CHUS_PT_GTVtot'); train2 = struct2cell(train2); train2 = train2{1};
test1 = load('AertsSign_HMR_PT_GTVtot'); test1 = struct2cell(test1); test1 = test1{1};
test2 = load('AertsSign_CHUM_PT_GTVtot'); test2 = struct2cell(test2); test2 = test2{1};
signTrain.PET = [train1.Data;train2.Data]; signTrain.PET = signTrain.PET(permTrain,:);
signTest.PET = [test1.Data;test2.Data]; signTest.PET = signTest.PET(permTest,:);
train1 = load('AertsSign_HGJ_CT_GTVtot'); train1 = struct2cell(train1); train1 = train1{1};
train2 = load('AertsSign_CHUS_CT_GTVtot'); train2 = struct2cell(train2); train2 = train2{1};
test1 = load('AertsSign_HMR_CT_GTVtot'); test1 = struct2cell(test1); test1 = test1{1};
test2 = load('AertsSign_CHUM_CT_GTVtot'); test2 = struct2cell(test2); test2 = test2{1};
signTrain.CT = [train1.Data;train2.Data]; signTrain.CT = signTrain.CT(permTrain,:);
signTest.CT = [test1.Data;test2.Data]; signTest.CT = signTest.CT(permTest,:);
for o = 1:nOutcomes
    training.(nameOutcomes{o}).sign = signTrain;
    testing.(nameOutcomes{o}).sign = signTest;
end
cd(pathCR), save('training','training'), save('testing','testing')
clear training testing


% FOR RF FOLDER
for o = 1:nOutcomes
    
    % FOR PET
    cd(pathLR), load('training'), load('testing'), pause(1)
    trainData = training.(nameOutcomes{o}).sign.PET; testData = testing.(nameOutcomes{o}).sign.PET;
    clear training testing, pause(1)
    cd(pathRF), load('training'), load('testing'), pause(1)
    PT_Energy = trainData(:,1); PT_Compactness = trainData(:,2); PT_GLN = trainData(:,3); PT_GLN_HLH = trainData(:,4);
    training.textures.([nameOutcomes{o},'Sign']).PET = table(PT_Energy,PT_Compactness,PT_GLN,PT_GLN_HLH);
    PT_Energy = testData(:,1); PT_Compactness = testData(:,2); PT_GLN = testData(:,3); PT_GLN_HLH = testData(:,4);
    testing.textures.([nameOutcomes{o},'Sign']).PET = table(PT_Energy,PT_Compactness,PT_GLN,PT_GLN_HLH);
    save('training','training'), save('testing','testing'), pause(2)
    clear training testing
    
    % FOR CT
    cd(pathLR), load('training'), load('testing'), pause(1)
    trainData = training.(nameOutcomes{o}).sign.CT; testData = testing.(nameOutcomes{o}).sign.CT;
    clear training testing, pause(1)
    cd(pathRF), load('training'), load('testing'), pause(1)
    CT_Energy = trainData(:,1); CT_Compactness = trainData(:,2); CT_GLN = trainData(:,3); CT_GLN_HLH = trainData(:,4);
    training.textures.([nameOutcomes{o},'Sign']).CT = table(CT_Energy,CT_Compactness,CT_GLN,CT_GLN_HLH);
    CT_Energy = testData(:,1); CT_Compactness = testData(:,2); CT_GLN = testData(:,3); CT_GLN_HLH = testData(:,4);
    testing.textures.([nameOutcomes{o},'Sign']).CT = table(CT_Energy,CT_Compactness,CT_GLN,CT_GLN_HLH);
    save('training','training'), save('testing','testing'), pause(2)
    clear training testing
    
end
cd(startpath)
end