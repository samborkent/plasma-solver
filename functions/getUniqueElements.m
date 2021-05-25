function [atomUC, nAtomUC] = getUniqueElements(uc, ucRatio)
%GETUNIQUEELEMENTS Get unique elements present in target and their
% normalized amount

% Get elements in unit cell
atomUCTemp = [uc.ELEMENTS];

% Initialize atom arrays
nAtomUCTemp     = [];  % Array containing amount of each atom in the target

% Get number of atoms weighted by ratio between unit cells
for iAtom = 1 : numel(uc)
    nAtomUCTemp = [nAtomUCTemp  (uc(iAtom).AMOUNT .* ucRatio(iAtom))];
end

% Normalize number of atoms array
nAtomUCTemp = nAtomUCTemp ./ sum(ucRatio);
nAtomUCTemp = nAtomUCTemp ./ min(nAtomUCTemp);

% Get sorted unique atom array
[elementNumber, elementIndex, ~] = unique([atomUCTemp.NUMBER], 'sorted');

% Check if there is oxygen (Z = 8)
if ~isempty(elementNumber(elementNumber == 8))
    % Get oxygen index
    oxygenIndex = elementIndex(elementNumber == 8);
    
    % Remove oxygen index
    elementIndex(elementNumber == 8) = [];
    
    % Add back oxygen index at the start
    elementIndex = [oxygenIndex; elementIndex];
end

% Insert unique atoms in array
atomUC = atomUCTemp(elementIndex');

% Initialize amount of atom array
nAtomUC = zeros(1, numel(elementNumber));

% Look through unique atoms
for i = 1 : numel(elementNumber)    
    % Insert the calculated ratio of each atom in atom amount array
    nAtomUC(i) = sum(nAtomUCTemp([atomUCTemp.NUMBER] == atomUC(i).NUMBER));
end

end

