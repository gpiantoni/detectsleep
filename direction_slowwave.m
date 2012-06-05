function [sw] = direction_slowwave(cfg, sw)
%DIRECTION_SLOWWAVE based on the delays, it shows the traveling of SW
% Use as:
%   [sw] = direction_slowwave(cfg, sw)
%
% CFG
%  .layout: layout of the electrodes
%
% (optional)
%  .dx: resolution of the grid (default 0.01)
%  .feedback: see ft_progress
%  .plot: show a plot for each traveling wave (default false)
%
% SW is the output of DETECT_SLOWWAVE
%
% This function adds a fields:
%   .streamline: the coordinates of the longest streamline
%   .origin: the starting point of the slow wave
%
% Part of DETECTSLEEP
% See also DETECT_SLOWWAVE DETECT_SPINDLE FIND_SPINDLES_PEAKFREQ

%---------------------------%
%-prepare input
%-----------------%
%-check cfg
if ~isfield(cfg, 'dx'); cfg.dx = 0.01; end
if ~isfield(cfg, 'plot'); cfg.plot = false; end
if ~isfield(cfg, 'feedback'); cfg.feedback = 'textbar'; end
%-----------------%

%-----------------%
%-layout
layout = cfg.layout;
%-----------------%
%---------------------------%

%---------------------------%
%-loop over slow waves
ft_progress('init', cfg.feedback)

for i = 1:numel(sw)

  ft_progress(i / numel(sw))
  
  %-----------------%
  %-location and timing of the negative peak
  [~, idelay, ilay] = intersect(sw(i).label, layout.label);
  
  eegx = layout.pos(ilay,1);
  eegy = layout.pos(ilay,2);
  eegz = sw(i).delay(idelay);
  %-----------------%
  
  %-----------------%
  %-create grid
  eegxmin = round(min(eegx)/cfg.dx) * cfg.dx;
  eegxmax = round(max(eegx)/cfg.dx) * cfg.dx;
  eegymin = round(min(eegy)/cfg.dx) * cfg.dx;
  eegymax = round(max(eegy)/cfg.dx) * cfg.dx;
  
  [grid_x, grid_y] = meshgrid(eegxmin:cfg.dx:eegxmax, eegymin:cfg.dx:eegymax);
  %-----------------%
  
  %-----------------%
  %-interpolation and gradient
  V = griddata(eegx, eegy, eegz, grid_x, grid_y);
  
  [grad1, grad2] = gradient(V);
  %-----------------%
  
  %-----------------%
  %-streamlines
  % we start a streamline at each point of the grid, it could be faster if
  % you skipped some points
  start_x = grid_x;
  start_y = grid_y;
  
  swstream = stream2(grid_x, grid_y, grad1, grad2, start_x, start_y);
  %-----------------%
  
  %-----------------%
  %-longest streamline
  lenstr = streamlength(swstream);
  
  [~, longest] = max(lenstr);
  trvl = swstream{longest};
  
  firstnan = find(isnan(trvl(:,1)), 1);
  trvl(firstnan:end,:) = [];
  
  sw(i).streamline = trvl;
  sw(i).origin = trvl(1,:);
  %-----------------%
  
  %-----------------%
  %-feedback
  if cfg.plot
    
    h = figure;
    subplot(1,2,1)
    ft_plot_lay(layout, 'label', 'no', 'point', 'no', 'box', 'no')
    hold on
    surf(grid_x, grid_y, V, 'EdgeColor', 'none');
    colorbar
    plot3(eegx, eegy, ones(size(eegx)), '.k') % you need plot3 bc surf is actually 3d and it would cover the dots
    
    %-------%
    %-get longest streamlines
    slenstr = sort(lenstr);
    toplot = lenstr >= slenstr(end-round(numel(slenstr) / 90)); % 90% percentile roughly
    %-------%
    
    %-------%
    %-plot it
    subplot(1,2,2)
    contour(grid_x, grid_y, V, 20)
    hold on
    quiver(grid_x, grid_y, grad1, grad2, 'k')
    h1 = streamline(swstream(toplot));
    set(h1,'color','r')
    axis equal
    %-------%
    
    waitfor(h)
  end
  %-----------------%
  
end

ft_progress('close')
%---------------------------%

%---------------------------------------------------------%
%-calculate the length of the SW streams------------------%
%---------------------------------------------------------%
function [lenstr] = streamlength(SWstr)

lenstr = zeros(numel(SWstr), 1);
for st = 1:numel(SWstr)
  firstnan = find(isnan(SWstr{st}(:,1)), 1);
  if isempty(firstnan); firstnan = 0; end
  lenstr(st) = firstnan;
end
%---------------------------------------------------------%