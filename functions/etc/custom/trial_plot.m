function trial_plot(data, including, threshold, savepath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trial_plot(data, including, threshold, savepath)
%
% Plot sensor timeseries in all trials vusialization of the
% inclusion/exclusion of trials. Original code by Simeon Wong.
%
% Inputs:
%   data                (structure) MEG data in FieldTrip format
%   including           (array)     Logical array indicating whether or not
%                                   a given trial should be included in the
%                                   analyses
%                                   (n x 1)
%                                       n   Number of trials
%   threshold           (float)     Threshold used to determine inclusion
%   savepath            (string)    Prefix for outputted figures
%
% Outputs:
%   "savepath*"         (files)     Figures visualizing trials
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify some options
p.Results.subplotsX = 5;
p.Results.subplotsY = 4;
p.Results.position = [10 10 1800 1000];

% Perpare plot
f = figure;
f.Position = p.Results.position;
f.PaperPositionMode = 'auto';

% Get the lobes for colouring purposes
lobeColors = lines(length(data.label));

% Determine number of plots and subplots
nSubplots = p.Results.subplotsX * p.Results.subplotsY;
nPlots = ceil(size(data.trial, 2) / nSubplots);

% Iterate over the number of plots
for pp = 1:nPlots
    
  figure(f)
  clf(f)
  
  % Iterate over the subplots
  for sp = 1:nSubplots
      
    % Break if we're done with trials
    tt = (pp-1) * nSubplots + sp;
    if tt > size(data.trial, 2)
      break
    end
    
    % Define the axis
    ax = subplot(p.Results.subplotsY, p.Results.subplotsX, sp);
    
    % Iterate over the lobes and plot
    for ll = 1:length(data.label)
      h = plot(data.time{1}, data.trial{tt}(ll, :), 'Color', lobeColors(ll,:));
      h.DataTipTemplate.DataTipRows(end+1) = data.label{ll};
      hold on
    end
    
    % Plot thresholds
    plot(data.time{1}([1, end]), [1 1] * threshold, ':k');
    plot(data.time{1}([1, end]), [1 1] * (-threshold), ':k');
    
    % Axis labels
    xlabel('Time (s)');
    ylabel('Amplitude (T)');
    title(['Trial ',num2str(tt)]);
    
    % Set the limits    
    ylim([-threshold,threshold] * 1.5);
    xlim([data.time{1}(1), data.time{1}(end)]);
    
    % Shade if including
    if including(tt) == 0
      ax = gca;
      ax.Color = [1.0 0.95 0.95];
    end
    f.InvertHardcopy = 'off';
    f.Color = 'white';
    
  end
  
  % Save
  saveas(f, [savepath '_' num2str(pp) '.png']);
  
end


% close(f);

end