%% Execute this script with fastAUC.cpp in the current directory
%Compile
mex fastAUC.cpp

%Test - needs the statistics toolbox
labels   = [1,1,1,1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,1,1,1,1];
scores   = rand(1,numel(labels));
posclass = 1;
AUC = fastAUC(labels,scores,posclass);
[~,~,~,AUC2] = perfcurve(labels,scores,posclass);

if(AUC == AUC2)
    fprintf('Test succeeded, AUC=%.4f\n',AUC);
else
    fprintf('Test failed! %.4f != %.4f\n',AUC,AUC2);
end