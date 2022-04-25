function massSpec = setupMassSpec(massSpecName)

switch massSpecName

    case "PhoenixKansas_1e12"
    massSpec.collectorWidthMM = 1;
    massSpec.effectiveRadiusMagnetMM = 540;
    massSpec.faradayNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
    massSpec.ionCounterNames = ["PM", "SEM"];
    massSpec.amplifierResistance = 1e12;

    case "PhoenixKansas_1e11"
    massSpec.collectorWidthMM = 1;
    massSpec.effectiveRadiusMagnetMM = 540;
    massSpec.faradayNames = ["L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
    massSpec.ionCounterNames = ["PM", "SEM"];
    massSpec.amplifierResistance = 1e12;

end