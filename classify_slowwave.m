function sw = classify_slowwave(cfg, sw)
%DIRECTION_SLOWWAVE based on the delays, it shows the traveling of SW
% Use as:
%   [sw] = classify_slowwave(cfg, sw)
%
% CFG
%  .fun:
%       'stepwise': it counts each step and adds them up
%       'beginend':
%
%  .mintrvl: min distance to travel (default 0)
%  .else: what to do with waves which don't pass mintrvl
%         ('back2front' 'front2back' 'none')
%
%-'beginend': check where the beginning or the end are more frontal
%
% SW is the output of DETECT_SLOWWAVE
%
% Part of DETECTSLEEP
% See also DETECT_SLOWWAVE DETECT_SPINDLE FIND_SPINDLES_PEAKFREQ

%---------------------------%
%-prepare input
%-----------------%
%-check cfg
if ~isfield(cfg, 'mintrvl'); cfg.mintrvl = 0; end
if ~isfield(cfg, 'else'); cfg.else = 'none'; end
if ~isfield(cfg, 'feedback'); cfg.feedback = 'textbar'; end
%-----------------%
%---------------------------%

%---------------------------%
%-loop over slow waves
ft_progress('init', cfg.feedback)

for i = 1:numel(sw)
  
  ft_progress(i / numel(sw))
  [f2b b2f] = feval(cfg.fun, sw(i).streamline);
  
  %-----------------%
  %-classify by counting steps
  if b2f - f2b > cfg.mintrvl
    sw(i).type = 'back2front';
    sw(i).param = abs(b2f - f2b);
    
  elseif f2b - b2f > cfg.mintrvl
    sw(i).type = 'front2back';
    sw(i).param = abs(b2f - f2b);
    
  else
    
    sw(i).type = cfg.else;
    
    switch cfg.else
      case 'back2front';
        sw(i).param = abs(b2f - f2b);
        
      case 'front2back';
        sw(i).param = abs(b2f - f2b);
        
      otherwise
        sw(i).param = 0;
        
    end
    
  end
  %-----------------%
  
end

ft_progress('close')
%---------------------------%

%---------------------------%
%-STEPWISE (if the longest stream goes towards the back)
%---------------------------%
function [f2b b2f] = stepwise(stream)

Ydiff = diff(stream(:,2));
f2b = numel(find(Ydiff < 0));
b2f = numel(find(Ydiff > 0));
%---------------------------%

%---------------------------%
%-BEGINEND (compare start and end point: absolute difference)
%---------------------------%
function [f2b, b2f] = beginend(stream)

f2b = stream(1,2);
b2f = stream(end,2);
%---------------------------%
