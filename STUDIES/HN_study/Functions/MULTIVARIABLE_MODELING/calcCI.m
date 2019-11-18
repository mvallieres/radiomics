function CI = calcCI(resp,time,censoring)
% resp: Hazards ratio for Cox regression
% Y: Time-to-event (the higher, the better)

% FLIP THE RESPONSE. If Resp1 > resp2 && Y1 > Y2, add + 1. IF Y1<Y2 and Y1
% is censored, discard.
% ________________________________________________________________________


% INITIALIZATION
nInst = numel(time);
resp = - resp; % Flipping the hazards for easier comparison


% GETTING ALL PAIRS OF RESPONSES
pairs = combnk(1:nInst,2);
resp = resp(pairs); 
time = time(pairs);
censoring = censoring(pairs);


% REMOVING UNUSABLE DATA

% 1. Both patients lived: unusable data
sumCens = censoring(:,1) + censoring(:,2);
out = find(sumCens == 2);
resp(out,:) = []; time(out,:) = []; censoring(out,:) = [];

% 2. One patient died, but the follow-up of the other is not long enough (did not outlive the first patient): unusable data
valCheck = (censoring(:,2) - censoring(:,1)) .* (time(:,2) - time(:,1));
out = find(valCheck < 0);
resp(out,:) = []; time(out,:) = []; censoring(out,:) = [];

% 3. Both patients died at the same time: unusable data
valCheck1 = time(:,2) - time(:,1); valCheck2 = censoring(:,2) - censoring(:,1); % (all 1,1 censored pairs have been removed by now, leaving only 0,1 or 0,0)
out = find(valCheck1 == 0 & valCheck2 == 0);
resp(out,:) = []; time(out,:) = []; censoring(out,:) = [];


% COMPUTING TIES IN RESPONSES
diffResp = resp(:,2) - resp(:,1);
indTies = find(diffResp == 0);
nTies = numel(indTies);
resp(indTies,:) = []; time(indTies,:) = []; censoring(indTies,:) = [];


% CALCULATING CI
valCheck = (time(:,2) - time(:,1)) .* (resp(:,2) - resp(:,1)); nCheck = numel(valCheck);
indGood = find(valCheck > 0); nGood = numel(indGood);
CI = (nGood + 0.5*nTies)/(nCheck + nTies);
if isnan(CI)
    CI = 0.5;
end

end