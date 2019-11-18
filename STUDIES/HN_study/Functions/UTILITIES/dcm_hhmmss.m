function [totSec hh mm ss fract] = dcm_hhmmss(dateStr)
if ~ischar(dateStr)
    dateStr = num2str(dateStr);
end
hh = str2double(dateStr(1,1:2));
mm = str2double(dateStr(1,3:4));
ss = str2double(dateStr(1,5:6));
fract_start = find(dateStr,'.');
fract = str2double(dateStr(1,fract_start+1:end));
totSec = hh*60*60 + mm*60 + ss;
end