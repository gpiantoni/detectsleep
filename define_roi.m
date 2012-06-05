function roi = define_roi(cfg, layout)
%DEFINE_ROI define ROI
% Use as:
%  roi = define_roi(cfg, layout)
% 
% CFG
%  .roi(1).name: string with the roi name
%  .roi(1).center: center with x-y coordinates
%  .roi(1).radius: distance of the electrodes (default: 0.2)
% The units are the same as for the layout.
%
%  .feedback: logical, to plot feedback of ROI
% 
% LAYOUT: output of ft_prepare_layout
%
% ROI: a struct with fields
%  .name: name of electrode to average
%  .chan: cell with labels
% ROI can then be used in detect_slowwave and detect_spindle
% 
% Part of DETECTSLEEP

%---------------------------%
%-default cfg
if ~isfield(cfg, 'feedback'); cfg.feedback = false; end
for i = 1:numel(cfg.roi)
  if ~isfield(cfg.roi(i), 'radius') || isempty(cfg.roi(i).radius)
    cfg.roi(i).radius = 0.2; 
  end
end
%---------------------------%

%---------------------------%
%-loop over ROI
roi = [];
for i = 1:numel(cfg.roi)
  roi(i).name = cfg.roi(i).name;
  
  x0 = cfg.roi(i).center(1);
  y0 = cfg.roi(i).center(2);
  
  inside = sqrt((layout.pos(:,1) - x0).^2 + (layout.pos(:,2) - y0).^2) <= cfg.roi(i).radius;
  roi(i).chan = layout.label(inside);

end
%---------------------------%

%---------------------------%
%-feedback
if cfg.feedback
  
  %-------%
  %-parameters
  circle = @(r) sqrt(r^2 - ((-r:r/100:r).^2));
  colors = {'b' 'r' 'k' 'g' 'c' 'm' 'y'};
  %-------%
  
  figure
  ft_plot_lay(layout, 'box', 'no', 'point', 'no')
  hold on
  
  %-----------------%
  %-loop over ROI
  for i = 1:numel(cfg.roi)
    
    i_c = mod(i - 1, numel(colors)) + 1; % index for colors
    
    %-------%
    %-plot good labels
    [~, ilay] = intersect(layout.label, roi(i).chan);
    plot(layout.pos(ilay,1), layout.pos(ilay,2), ['.' colors{i_c}])
    %-------%
    
    %-------%
    %-plot circle
    r = cfg.roi(i).radius;
    x0 = cfg.roi(i).center(1);
    y0 = cfg.roi(i).center(2);
  
    uppercircle = y0 + circle(r);
    lowercircle = y0 - circle(r);
    
    plot((-r:r/100:r) + x0, uppercircle, colors{i_c})
    plot((-r:r/100:r) + x0, lowercircle, colors{i_c})
    %-------%
    
  end
  %-----------------%
  
end
%---------------------------%

