function setup_paths()
%SETUP_PATHS Add the source folder and CasADi to the MATLAB path.
%   CasADi (3.8.2 was used for the paper) is taken from the MATLAB path if it
%   is already there, otherwise from the folder casadi/ next to this file.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'));
if isempty(which('casadi.SX'))
    local = fullfile(root,'casadi');
    if isfolder(local), addpath(local); end
end
if isempty(which('casadi.SX'))
    error('setup_paths:CasADi', ['CasADi was not found. Download the MATLAB release of CasADi ' ...
        '(https://web.casadi.org/get/) and unzip it into %s, or add it to the MATLAB path.'], ...
        fullfile(root,'casadi'));
end
end
