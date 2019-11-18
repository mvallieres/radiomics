function calcAllFeatureImportanceRF_SUPP(pathRF,nPerms,seed,matlabPATH)

% VARYING PARAMETERS, DEPENDENT ON FINAL CHOICE OF BEST RESULTS (one Rad/Vol set for each outcome + clinical variables)
fSetNames = {'PETCTclinic','CTclinic','clinic'}; % For Locoregional, Distant and Death outcome, respectively
nameOutcomes = {'Locoregional','Distant','Death'};
nOutcomes = numel(nameOutcomes);


% INITIALIZATION
startpath = pwd;
time = 60;
cd(pathRF), load('training')
mkdir('batchLog_RFperm'), cd('batchLog_RFperm'), pathBatch = pwd;


% PRODUCE BATCH COMPUTATIONS
batch = 0;
save('workspace','pathRF','nPerms','nameOutcomes','fSetNames','seed')
for o = 1:nOutcomes
    batch = batch + 1;
    nameScript = ['batch',num2str(batch),'_script.m'];
    fid = fopen(nameScript,'w');
    fprintf(fid,'load(''workspace'')\n');
    fprintf(fid,['[percentAUCchange,varNames,pVal] = calcFeatureImportanceRF_HN(pathRF,nPerms,nameOutcomes{',num2str(o),'},fSetNames{',num2str(o),'},seed);\n']);
    fprintf(fid,['save(''result_batch',num2str(batch),''',''percentAUCchange'',''varNames'',''pVal'')\n']);
    fprintf(fid,['system(''touch batch',num2str(batch),'_end'');\n']);
    fprintf(fid,'clear all');
    fclose(fid);
    system([matlabPATH,' -nojvm -nodisplay -nodesktop -nosplash < ',nameScript,' >& ',nameScript(1:end-1),'log &']);
end

% WAITING LOOP
waitBatch(pathBatch,time,nOutcomes)
delete('workspace.mat')


% GROUP RESULTS IN THE 'training' STRUCT
cd(pathBatch)
batch = 0;
variableImportance = struct;
for o = 1:nOutcomes
    batch = batch + 1;
    load(['result_batch',num2str(batch)]) % Variables 'percentAUCchange', 'varNames' and 'pVal' gets out of there
    variableImportance.(nameOutcomes{o}).(fSetNames{o}).percentAUCchange = percentAUCchange;
    variableImportance.(nameOutcomes{o}).(fSetNames{o}).varNames = varNames;
    variableImportance.(nameOutcomes{o}).(fSetNames{o}).pVal = pVal;
    delete(['result_batch',num2str(batch),'.mat'])
    clear percentAUCchange varNames pVal
end
cd(pathRF), save('testingVariableImportance','variableImportance')

cd(startpath)
end