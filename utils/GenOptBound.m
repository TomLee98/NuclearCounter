function [lb, ub, nuclear] = GenOptBound(animal_info)
marker = upper(string(animal_info.marker));
driver = upper(string(animal_info.driver));
stage = upper(string(animal_info.stage));

switch marker
    case "MCHERRY"
        switch driver
            case "OK107"
                %     req   ssy   sim   avg   gam
                lb = [0.85; 0.90; 0.70; 0.00; 0.00];
                ub = [0.99; 0.99; 1.00; 0.10; 0.75];
                switch stage
                    case "L1"
                        nuclear = struct('radius', 2.1, ...
                                         'sigma', 0.25, ...
                                         'color',"r");
                    case "L2"
                        nuclear = struct('radius', 2.1, ...
                                         'sigma', 0.25, ...
                                         'color',"r");
                    case "L3"
                        throw("Unregistered Drosophila develop stage.");
                    case "ADULT"
                        throw("Unregistered Drosophila develop stage.");
                    otherwise
                        throw("Unknown Drosophila develop stage.");
                end
            case "ORCO"
                % TODO: parameters range modified after supervised learning
                %     req   ssy   sim   avg   gam
                lb = [0.85; 0.90; 0.50; 0.02; 0.00];
                ub = [0.99; 0.99; 1.00; 0.10; 0.95];
                switch stage
                    case "L1"
                        nuclear = struct('radius', 2.2, ...
                                         'sigma', 0.3, ...
                                         'color',"r");
                    case "L2"
                        throw("Unregistered Drosophila develop stage.");
                    case "L3"
                        throw("Unregistered Drosophila develop stage.");
                    case "ADULT"
                        throw("Unregistered Drosophila develop stage.");
                    otherwise
                        throw("Unknown Drosophila develop stage.");
                end
            otherwise
                throw("Unknown driver genes on Drosophila.");
        end
    otherwise
        throw("Unknown marker genes on Drosophila.");
end
end
