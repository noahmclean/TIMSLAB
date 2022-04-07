function [G, modelMasses] = assembleG(magnetMasses, massSpec, modelMassRange)

minModelMass = modelMassRange(1);
maxModelMass = modelMassRange(2);

averageMassAMU = mean(magnetMasses);
collectorWidthAMU = averageMassAMU / massSpec.effectiveRadiusMagnetMM ...
                    * massSpec.collectorWidthMM;

nMagnetMasses = length(magnetMasses);

% mass range that is inside the collector at any measured magnetMass
collectorLimits = magnetMasses + [-collectorWidthAMU, collectorWidthAMU]/2;

deltaMagnetMass = mean(magnetMasses(2:end-1)-magnetMasses(1:end-2));
nModelMasses = ceil(0.95*(maxModelMass-minModelMass)/deltaMagnetMass);
modelMasses = linspace(minModelMass, maxModelMass, nModelMasses);
deltaModelMass = modelMasses(2)-modelMasses(1);

% make a new convolution matrix, now with trapezoidal rule  
G = zeros(nMagnetMasses, nModelMasses);
for iMass = 1:nMagnetMasses % a row for each manget mass

    % massesInCollector are *model* masses
    massesInCollector = collectorLimits(iMass,1) <= modelMasses & ...
                        modelMasses <= collectorLimits(iMass,2);
    
    firstMassIndexInside = find(massesInCollector,1,'first');
    lastMassIndexInside  = find(massesInCollector,1,'last');
    G(iMass, firstMassIndexInside + 1: lastMassIndexInside -1) = deltaModelMass;
    G(iMass, [firstMassIndexInside, lastMassIndexInside]) = deltaModelMass/2;

    % add in fractional segments at edges of peak/collector overlap
    if firstMassIndexInside > 1 % if not on front edge
        
        m1 = modelMasses(firstMassIndexInside-1);
        m2 = modelMasses(firstMassIndexInside);
        mb = collectorLimits(iMass,1);
        G(iMass, firstMassIndexInside-1) = (m2-mb)^2/(2*(m2-m1));
        G(iMass, firstMassIndexInside) = G(iMass, firstMassIndexInside) + ...
                                    (m2-mb)*(m2-2*m1+mb)/(2*(m2-m1));

    end % if firstMassIndexInside > 1 (not on front edge)

    if lastMassIndexInside < nModelMasses % if not on back edge
        
        m1 = modelMasses(lastMassIndexInside);
        m2 = modelMasses(lastMassIndexInside+1);
        mb = collectorLimits(iMass,2);
        G(iMass, lastMassIndexInside) = G(iMass, lastMassIndexInside) + ...
                                (m1 - mb)*(m1-2*m2+mb)/(2*(m2-m1));
        G(iMass, lastMassIndexInside+1) = (mb-m1)^2/(2*(m2-m1));

    end % if firstMassIndexInside > 1 (not on back edge)

end

% note: this is a nicer/more accurate trapezoidal rule G
