classdef pt
    %PT Periodic Table data
    %   Contains atomic masses
    
    properties
        atomicMass double
    end
    
    methods
        function e = pt(isotopeName)
            %PT Construct an instance of this class
            %   with atomic masses
            e.atomicMass = isotopeName;
        end

        function m = mass(isotopeName)
            m = isotopeName.atomicMass;
        end

    end

    enumeration
        Sr84 (83.913425)
        Sr86 (85.9092607309)
        Sr87 (86.9088774970)
        Sr88 (87.9056122571)
        Rb85 (84.911789738)
        Rb87 (87.91131559)
    end

end

