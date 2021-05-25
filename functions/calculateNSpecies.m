function nSpecies = calculateNSpecies(atomUC, bg, nOxidesPerElement)
%CALCULATENSPECIES Calculate the number of possible species

% If there is oxygen in the target and in the background gas
if sum([bg.NUMBER] == 8) && sum([atomUC.NUMBER] == 8)
    nSpecies = 1 + numel(atomUC) + nOxidesPerElement*(numel(atomUC) - 1);
% If there is only oxygen in the background gas
elseif sum([bg.NUMBER] == 8)
    nSpecies = 1 + numel(atomUC) + nOxidesPerElement*nElements;
% If no oxygen is present in the background
else
    nSpecies = 1 + numel(atomUC);
end

end

