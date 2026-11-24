function W = beamform(sourcemodel, headmodel, LF, Cinv, method, nai)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W = beamform(cfg)
%
% Function to perform beamforming for the CRYO or OPM data in "opm format".
% Original code was provided by Dr. Ryan Hill.
%
% Inputs:
%	sourcemodel (structure) Fieldtrip sourcemodel
%	headmodel   (structure) Fieldtrip headmodel
%	LF          (array)     Lead fields from cryo/opm_create_leadfield
%                           (n x 3 x m)
%                            	n   Number of sensors
%                               m   Number of sources
%	Cinv        (array)     Inverse of the covariance matrix
%                           (n x n)
%   method      (string)    'both': Both orientations are kept
%                           'optimal': Optimal combination to maximize SNR
%                           (Sekihara et al., 2004)
%
% Outputs:
%   W           (array)     Beamformer weights
%                           (n x m)
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the headmodel and convert to meters to be consistent across
% leadfield and beamforming functions
sourcemodel_m = ft_convert_units(sourcemodel, 'm');
headmodel_m = ft_convert_units(headmodel, 'm');

% Extract the "origin" (mean position) of the headmodel (in m)
origin = mean(headmodel_m.bnd.pos);

% Check method to initialize the size of the weights
if strcmp(method, 'both')
    W = zeros(size(LF, 1), size(sourcemodel_m.pos, 1), 2);
elseif strcmp(method, 'optimal')
    W = zeros(size(LF, 1), size(sourcemodel_m.pos, 1), 1);
else
    error('Incorrect specification of ''method'': must be ''both'' or ''optimal''');
end

% Iterate over the source positions and get the beamformer weights
for ii = 1:size(sourcemodel_m.pos, 1)
    
    % Find vector to the point
    
    % Project to two tangential fields
    R = sourcemodel_m.pos(ii, :) - origin;
    [etheta,ephi] = calctangent(R);
    tan_ors = [etheta;ephi];
    ltan = LF(:,:,ii)*tan_ors';
    
    % Check method
    if strcmp(method, 'both')
        
        % Use both tangential fields
        lopt = ltan;
        
        % 
        
    else
        
        % Check if using NAI
        if nai
            
            % Optimal orientation calculation for unit-noise gain
            % beamformer
            
            % Singular value decomposition
            [v, d] = svd(pinv(ltan'*Cinv*ltan)*(ltan'*Cinv^2*ltan));
                        
        else
        
            % Get the optimal combination to maximize the SNR - Sekihara et al.
            % 2004

            % Singular value decomposition
            [v, d] = svd(ltan'*Cinv*ltan);
            
        end

        % Find the minimum
        [~, id] = min(diag(d));

        % Select this orientation and turn to a nAm amplitude
        lopt = ltan*v(:, id);
            
    end
    
    % Now, compute the weight
    for jj = 1:size(lopt, 2)
        if nai
            W(:, ii, jj) = Cinv*lopt(:, jj)/sqrt(lopt(:, jj)'*Cinv^2*lopt(:, jj));
        else
            W(:, ii, jj) = Cinv*lopt(:, jj)/(lopt(:, jj)'*Cinv*lopt(:, jj));
        end
    end
    
end


