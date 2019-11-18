function [cellIndex] = findUIDcellIndex(uidString,cellString)
% WRITTEN BY MARTIN: CREATE HEADER!

cellIndex = find(strcmp(uidString,cellString));
if isempty(cellIndex)
    cellIndex = numel(cellString) + 1; % Not present in cellString, so create a new position in the cell for the new UID
end

end