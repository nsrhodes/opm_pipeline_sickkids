function parsave(fname, x, str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parsave(fname, x, str)
%
% Function to save a variable within a parpool loop
%
% Inputs:
%   x               (variable)  Variable to save
%   str             (string)    Name to save the variable under
%
% Outputs:
%   "fname"         (file)      Saved file
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    eval([str '= x;']);
    save(fname, '-v7.3', str)
    
end