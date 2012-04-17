function [SP] = detect_spindle(cfg, data)
%DETECT_SPINDLE detect spindles using simple hilbert
% based on 2 criteria:
% 1. bandpass hilbert is more than thr
% 2. duration above thr is longer than dur
%
% Use as:
%    [SP] = detect_spindle(cfg, data)
%
% cfg
%  .roi(1).name = name of electrode to average
%  .roi(1).chan = chan index or cell with labels (better)
%
%  (optional)
%  .filter = [11 15]; (can be empty)
%  .zeropad = 1; (if filter, zero padding to avoid edge artifacts)
%  .thr = 3; (threshold for bandpassed hilbert abs)
%  .dur = [.5 1]; (duration above the threshold)
%  if cfg.dur = []; then you can specify cfg.thrB, the second threshold like in Ferrarelli et al. 2007 but it's not adaptive
%  .feedback = 'textbar' (see ft_progress)
%
% data
%    data in fieldtrip format, see also sleep2ft
%    can be a filename, it'll read variable 'data'
%
% SP
%    struct for each detected spindle
%  .trl = trial it belongs to
%
%  .begsp_itrl = first point above threshold in samples, from beginning of the trial
%  .begsp_iabs = first point above threshold in samples, based on sampleinfo
%  .begsp_time = first point above threshold in seconds, based on data.time
%
%  .max        = highest value
%  .maxsp_itrl = highest value in samples, from beginning of the trial
%  .maxsp_iabs = highest value in samples, based on sampleinfo
%  .maxsp_time = highest value in seconds, based on data.time
%
%  .endsp_itrl = last point above threshold in samples, from beginning of the trial
%  .endsp_iabs = last point above threshold in samples, based on sampleinfo
%  .endsp_time = last point above threshold in seconds, based on data.time
%
%  .energytot  = total amount of energy (normalized by sampling rate: min is cfg.dur(1) * cfg.thr)
%  .energysec  = total amount of energy divided by time (in seconds)
%

%---------------------------%
%-prepare input
%-----------------%
%-check cfg
%-------%
%-defaults
if ~isfield(cfg, 'filter'); cfg.filter = [11 15]; end
if ~isfield(cfg, 'zeropad'); cfg.zeropad = 1; end
if ~isfield(cfg, 'dur'); cfg.dur = [.5 1]; end
if ~isfield(cfg, 'feedback'); cfg.feedback = 'textbar'; end

if ~isempty(cfg.dur) % by default use duration method
  if ~isfield(cfg, 'thr'); cfg.thr = 3; end
else
  if ~isfield(cfg, 'thr'); cfg.thr = 7; end % detection threshold needs to be higher
  if ~isfield(cfg, 'thrB'); cfg.thrB = 2; end % threshold for duration
end  
%-------%
%-----------------%

%-----------------%
%-check data
if ischar(data)
  load(data, 'data')
end
%-----------------%
%---------------------------%

%---------------------------%
%-prepare the data
%-----------------%
%-filtering with zero-padding
if ~isempty(cfg.filter)
  
  %-------%
  %-zero padding
  nchan = numel(data.label);
  pad = cfg.zeropad * data.fsample;
  datatime = data.time; % keep it for later
  
  for i = 1:numel(data.trial)
    data.trial{i} = [zeros(nchan,pad) data.trial{i} zeros(nchan,pad)];
    data.time{i} = (1:size(data.trial{i},2)) / data.fsample;
  end
  %-------%
  
  %-------%
  %-filter
  cfg1 = [];
  cfg1.hpfilter = 'yes';
  cfg1.hpfreq = cfg.filter(1);
  cfg1.hpfiltord = 4;
  
  cfg1.lpfilter = 'yes';
  cfg1.lpfreq = cfg.filter(2);
  
  cfg1.feedback = cfg.feedback;
  cfg1.hilbert = 'abs';
  
  data = ft_preprocessing(cfg1, data);
  %-------%
  
  %-------%
  %-remove padding
  for i = 1:numel(data.trial)
    data.trial{i} = data.trial{i}(:, (pad+1) : (end-pad));
  end
  
  data.time = datatime;
  %-------%
  
end
%-----------------%

%-----------------%
%-define roi
datroi = rmfield(data, {'trial' 'label'});
for r = 1:numel(cfg.roi)
  datroi.label{r,1} = cfg.roi(r).name;
  
  if iscell(cfg.roi(r).chan)
    [~, roichan] = intersect(data.label, cfg.roi(r).chan);
  else
    roichan = cfg.roi(r).chan;
  end
  %-------%
  %-average over channels
  for i = 1:numel(data.trial)
    datroi.trial{i}(r,:) = mean(data.trial{i}(roichan, :), 1);
  end
  %-------%
  
end
%-----------------%
%---------------------------%

%---------------------------%
%-detect spindle
spcnt = 0;
SP = [];

for r = 1:numel(datroi.label)
  
  for trl = 1:numel(datroi.trial);
    
    ft_progress(trl / numel(datroi.trial))
    x  = datroi.trial{trl}(r,:);
    
    %-----------------%
    %- 1. above threshold
    absp = find(x >= cfg.thr);
    if numel(absp) ==0
      continue
    end
    %-----------------%
    
    begsp = [1 find([1 diff(absp)] ~= 1)];
    endsp = [find([1 diff(absp)] ~= 1)-1 numel(absp)];
    
    for s = 1:numel(begsp)
      
      %-----------------%
      %- 2. second criterion
      if ~isempty(cfg.dur)
        %-------%
        %-duration above threshold
        dursp = (endsp(s) - begsp(s)) / data.fsample;
        isspindle = dursp >= cfg.dur(1) && dursp <= cfg.dur(2);
        
        begsp_itrl = absp(begsp(s));
        endsp_itrl = absp(endsp(s));
        %-------%
        
      else
        
        %-------%
        %-second threshold a la Ferrarelli 2007
        begsp_itrl = find( x(1: absp(begsp(s))) < cfg.thrB, 1, 'last');
        
        if ~isempty(begsp_itrl) 
          endsp_itrl = begsp_itrl + find( x((begsp_itrl+1):end) < cfg.thrB, 1, 'first');
          
          if ~isempty(endsp_itrl)
            isspindle = true;
          else
            isspindle = false;
          end
        else
          isspindle = false;
        end
        %-------%
      end
      
      if isspindle
        spcnt = spcnt + 1;
        
        [SP(spcnt).max, imax] = max(x(begsp_itrl:endsp_itrl)); % max point of the spindle
        
        SP(spcnt).energytot = sum(x(begsp_itrl:endsp_itrl)) / data.fsample;
        SP(spcnt).energysec = SP(spcnt).energytot / ((endsp_itrl - begsp_itrl) / data.fsample);
        SP(spcnt).trl = trl;
        SP(spcnt).begsp_itrl = begsp_itrl;
        SP(spcnt).begsp_iabs = SP(spcnt).begsp_itrl + datroi.sampleinfo(trl, 1);
        SP(spcnt).begsp_time = datroi.time{trl}(SP(spcnt).begsp_itrl);
        
        SP(spcnt).maxsp_itrl = begsp_itrl + imax -1;
        SP(spcnt).maxsp_iabs = SP(spcnt).maxsp_itrl + datroi.sampleinfo(trl, 1);
        SP(spcnt).maxsp_time = datroi.time{trl}(SP(spcnt).maxsp_itrl);
        
        SP(spcnt).endsp_itrl = endsp_itrl;
        SP(spcnt).endsp_iabs = SP(spcnt).endsp_itrl + datroi.sampleinfo(trl, 1);
        SP(spcnt).endsp_time = datroi.time{trl}(SP(spcnt).endsp_itrl);
        
      end
      %-----------------%
      
    end
    
  end
end

if ~isempty(SP)
  %-----------------%
  %-remove duplicates
  %-------%
  %-pure duplicates
  % (because of padding the same data is used in two consecutive trials, for example)
  % check if the values of two negpeak in absolute samples are the same
  dupl = [true diff([SP.maxsp_iabs]) ~= 0];
  SP = SP(dupl);
  %-------%
  
  %-------%
  %-similar, from different ROIs
  [~, SPi] = sort([SP.maxsp_iabs]);
  SP = SP(SPi);
  
  % don't use only negpeak
  dupl = [true diff([SP.maxsp_iabs]) > 50]; % TODO: specify this number
  SP = SP(dupl);
  %-------%
  %-----------------%
end
%---------------------------%

