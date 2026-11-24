function LF = opm_create_leadfield(grad, headmodel, sourcemodel, M, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LF = opm_create_leadfield(cfg)
%
% Function to create lead fields for the OPM. Modified code provided by Dr.
% Ryan Hill
%
% Inputs:
%	grad        (structure)     Fieldtrip gradiometer structure
% 	headmodel   (structure)     Fieldtrip headmodel
% 	sourcemodel (structure)     Fieldtrip sourcemodel
%   M           (array)         Mean field correction matrix
%                               (n x n)
%                                   n   Number of sensors
%
% Outputs:
%	LF          (array)         Lead fields 
%                               (n x 3 x m)
%                                   m   Number of sources
%                                   3   (x, y, z)
%       
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the headmodel, source model, and gradiometers, and convert to
% meters to be consistent
headmodel_m = ft_convert_units(headmodel, 'm');
sourcemodel_m = ft_convert_units(sourcemodel, 'm');
grad_m = ft_convert_units(grad, 'm');

% Create the leadfield
cfg_lf = [];
cfg_lf.grad = grad_m;
cfg_lf.headmodel = headmodel_m;
cfg_lf.sourcemodel = sourcemodel_m;
cfg_lf.units = 'm';
Lead_fields = ft_prepare_leadfield(cfg_lf);

% Reshape
LF = zeros(length(Lead_fields.label), 3, size(sourcemodel_m.pos, 1));
for n = 1:size(sourcemodel_m.pos, 1)
    
    % Assign leadfield if inside the head
    if Lead_fields.inside(n)
        LF(:, :, n) = Lead_fields.leadfield{n}; 
    else
        LF(:, :, n) = zeros(length(Lead_fields.label), 3);
    end
    
    % Account for mean field correction - note, if tSSS was used, 
    % this is autmatically set to the identity matrix so no effect.
    if exist('M', 'var')
        if strcmp(method, 'tAMM')
            for n = 1:size(LF,3)
                LF(:,:,n) = (M.Rout*M.Sin*pinv(M.Rout*M.Sin))*LF(:,:,n);
            end
        else
            LF(:, :, n) = M*LF(:, :, n);
        end
    end

        % for n = 1:size(LF,3)
        %     LF(:,:,n) = (M.Rout*S.Sin*pinv(M.Rout*S.Sin))*LF(:,:,n);
        % end
    
end
    

% Convert to nano-meters
LF = LF.*1e-9;
