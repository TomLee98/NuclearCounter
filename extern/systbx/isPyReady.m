function tf = isPyReady()
%ISPYREADY This function returns python environment status
% check the environment
MAT_PY_VALID_MAPPING = struct("R2020a",["2.7","3.6","3.7"],...
    "R2020b",["2.7","3.6","3.7","3.8"],...
    "R2021a",["2.7","3.7","3.8"],...
    "R2021b",["2.7","3.7","3.8","3.9"],...
    "R2022a",["2.7","3.7","3.8","3.9"],...
    "R2022b",["2.7","3.8","3.9","3.10"],...
    "R2023a",["3.8","3.9","3.10"], ...
    "R2023b",["3.9","3.10","3.11"]);

% get matlab version
mat_ver = string(version).extractBetween("(",")");

pe = pyenv;

if (pe.Version ~= "")
    if ismember(pe.Version, MAT_PY_VALID_MAPPING.(mat_ver))
        if ispc()
            % redirect output to avoid console display
            s_np = system("pip list | findstr numpy > ~.sftxt");
            s_tf  = system("pip list | findstr tifffile >> ~.sftxt");
        elseif isunix()
            s_np = system("pip list | grep numpy > ~.sftxt");
            s_tf  = system("pip list | grep tifffile >> ~.sftxt");
        else
            throw(MException("isPyReady:invalidOperationSystem", ...
                "Sorry, your operation system has not been supported."));
        end
        delete("~.sftxt");    % remove the redirect temporary file
        if (s_np==0) && (s_tf==0)
            % there are python and tifffile package, use tifffile replace
            % bio-formats java package
            tf = true;
        else
            tf = false;
            warning("savetiff:noTFOrNP", ...
                "Fast saving disabled because of numpy or tifffile package(s) lost.");
        end
    else
        tf = false;
        warning("savetiff:unmatchedPythonVersion","Fast saving disabled because of " + ...
            "python version is not match.");
    end
else
    tf = false;
    warning("savetiff:noPythonOrNotPC","Fast saving disabled because of " + ...
        "unsupported platform or python environment is not installed.");
end

end

