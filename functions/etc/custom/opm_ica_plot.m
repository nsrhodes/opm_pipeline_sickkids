function opm_ica_plot(comp, grad, save_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cryo_ica_plot(comp, grad, save_name)
%
% Function to visualize ICA components derived from CRYO data processed 
% using FieldTrip
%
% Inputs:
%   comp            (structure) FieldTrip structure containing the ICA
%                               components
%   grad            (structure) Fieldtrip gradiometer structure
%   save_name       (string)    Prefix for output figures
%
% Outputs:
%   "save_name*"    (files)     Figures visualizing ICA components
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of components
num_comp = size(comp.trial{1}, 1);

% Split into groups of 10
num_comp_per_plot = ceil(num_comp / ceil(num_comp / 10));
num_plots = ceil(num_comp / num_comp_per_plot);

% Make each plot
for pp = 1:num_plots
    
    % Setup figure
    hF = figure;
    hF.Units = 'pixels';
    hF.Position = [50 50 1200 900];
    padding = 10;
    rowheight = (hF.Position(4) - padding * 2) / num_comp_per_plot;
    topo_pos = @(n) round([padding, rowheight * (n-1) + padding*2, rowheight - padding * 2, rowheight - padding * 2]);
    ts_pos = @(n) round([rowheight + padding * 2, rowheight * (n-1) + padding*2, hF.Position(3) - rowheight - padding * 4, rowheight - padding * 2]);

    % Get the 95th percentile of all the data
    temp = zeros(size(comp.trial{1}, 1), size(comp.trial{1}, 2), length(comp.trial));
    for ii = 1:size(temp, 3)
        temp(:, :, ii) = comp.trial{ii};
    end
    m1 = prctile(temp(:), 95);
    m2 = prctile(temp(:), 5);
    m = max([abs(m1), abs(m2)]);
    
    for cc = 1:num_comp_per_plot
        kk = (pp-1)*num_comp_per_plot + cc;
        if kk > num_comp
            break;
        end
        
        % Make the topographical plot
        ax = axes();
        ax.Units = 'pixels';
        ax.Position = topo_pos(num_comp_per_plot+1-cc);
        ft_topoplotIC(struct('figure', hF, 'component', kk, 'grad', grad, 'comment', 'no', 'markersymbol', '.', 'title', ['C ', num2str(kk)]), comp);
        
        % Make the timeseries plot
        ax = axes();
        ax.Units = 'pixels';
        ax.Position = ts_pos(num_comp_per_plot+1-cc);
        hold on
        for tt = 1:min(length(comp.trial), 30)
            plot(comp.time{1}, comp.trial{tt}(kk,:));
            ylim(5*[-m, m]);
            xlim tight
        end
        
        annotation('textbox', 'String', ['C',num2str(kk)], 'LineStyle', 'none', 'Units', 'pixels', ...
            'Position', ts_pos(num_comp_per_plot+1-cc), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
            'FontWeight', 'bold', 'FontSize', 9, 'Margin', 2);
    end
    
    % save the plots
    saveas(hF, save_name(pp));
    
end

