function [SW] = detect_slowwave(cfg, data)
%DETECT_SLOWWAVE detect slow waves following Massimini, 2004
% based on 3 criteria:
% 1. less than negthr
% 2. distance between zero-crossing
% 3. peak-to-peak (if cfg.p2p is not empty)
% 
% Use as:
%    [SW] = detect_slowwave(cfg, data)
% 
% cfg
%  .roi(1).name: name of electrode to average
%  AND
%  .roi(1).chan: chan index or cell with labels (better)
%  roi can be defined using DEFINE_ROI
%
%  (optional)
%  .preproc = struct which is passed to ft_preprocessing ft_preprocessing(cfg.preproc, data)
%             If not specified, it's low-pass filter at 4 Hz and high-pass filter at 0.25 HZ
%             If empty, no preprocessing
%             See ft_preprocessing
%  .zeropad = 1; (if filter, zero padding to avoid edge artifacts)
%  .negthr = -40; (min uV to detect a slow wave)
%  .zcr = [.2 1]; (min and max distance between zero-crossing
%  .p2p = 75; (min uV peak-to-peak amplitude; if [], doesn't run this part)
%  .postzcr = 1; (max distance after zerocrossing, to be test for p2p)
%
%  (options for traveling)
%  .trvl.negthr: -30 (min uV to consider a channel having a slow wave,
%                 more liberal than for sw_detect)
%  .trvl.wnw: window around negative peak to look for each channel's
%             negative peak (two numbers in ms, [-.2 .2])
%  .trvl.dist: distance in s between the negative peaks of various
%              electrodes to consider them in the main wave
% 
%  .feedback = 'textbar' (see ft_progress)
%
% data
%    data in fieldtrip format, see also sleep2ft
%    can be a filename, it'll read variable 'data'
%
% SW
%    struct for each detected slow wave
%  .trl = trial it belongs to
%
%  .negpeak_val = value of the negative peak (uV)
%  .negpeak_itrl = negative peak in samples, from beginning of the trial
%  .negpeak_iabs = negative peak in samples, based on sampleinfo
%  .negpeak_time = negative peak in seconds, based on data.time
%
%  .begsw_itrl = 0-crossing down in samples, from beginning of the trial
%  .begsw_iabs = 0-crossing down in samples, based on sampleinfo
%  .begsw_time = 0-crossing down in seconds, based on data.time
%
%  .zcr_itrl = 0-crossing up in samples, from beginning of the trial
%  .zcr_iabs = 0-crossing up in samples, based on sampleinfo
%  .zcr_time = 0-crossing up in seconds, based on data.time
%
%  .pospeak_val = value of the positive peak (uV)
%  .pospeak_itrl = positive peak in samples, from beginning of the trial
%  .pospeak_iabs = positive peak in samples, based on sampleinfo
%  .pospeak_time = positive peak in seconds, based on data.time
%   (then cfg.p2p, if not empty, tests that pospeak_val - negpeak_val > cfg.p2p)
%
%  .channel = order of channels, based on their negative peaks
%  .label = same as channel, with labels
%  .delay = delay from first negative peak
%  .maxdelay = maxdelay
%
% Note that FASST uses a 90% criterion on the max slope. It's not reported
% in the article and I prefer not to include it. The results are therefore
% different
%
% Part of DETECTSLEEP
% See also DETECT_SLOWWAVE DETECT_SPINDLE FIND_SPINDLES_PEAKFREQ

%---------------------------%
%-prepare input
%-----------------%
%-check cfg
%-------%
%-defaults
if ~isfield(cfg, 'preproc');
  cfg.preproc.lpfilter = 'yes';
  cfg.preproc.lpfreq = 4;
  cfg.preproc.lpfiltord = 4;
  
  cfg.preproc.hpfilter = 'yes';
  cfg.preproc.hpfreq = 0.25;
  cfg.preproc.hpfiltord = 4;
  
end

%-------%
%-other parameters
if ~isfield(cfg, 'zeropad'); cfg.zeropad = 1; end
if ~isfield(cfg, 'negthr'); cfg.negthr = -40; end
if ~isfield(cfg, 'zcr'); cfg.zcr = [.2 1]; end
if ~isfield(cfg, 'p2p'); cfg.p2p = 75; end
if ~isfield(cfg, 'postzcr'); cfg.postzcr = 1; end
if ~isfield(cfg, 'feedback'); cfg.feedback = 'textbar'; end
if ~isfield(cfg, 'trvl'); cfg.trvl = []; end
if ~isfield(cfg.trvl, 'negthr'); cfg.trvl.negthr = -30; end
if ~isfield(cfg.trvl, 'wnw'); cfg.trvl.wnw = [-.05 .2]; end
if ~isfield(cfg.trvl, 'dist'); cfg.trvl.dist = 0.05; end
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
if ~isempty(cfg.preproc)
  
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
  cfg1 = cfg.preproc;
  cfg1.feedback = cfg.feedback;
  
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
%-define roip
mont = prepare_montage(cfg.roi, data.label);
datroi = ft_apply_montage(data, mont, 'feedback', 'none');
%-----------------%
%---------------------------%

%---------------------------%
%-detect slow wave
swcnt = 0;
SW = [];

for r = 1:numel(datroi.label)
  
  ft_progress('init', cfg.feedback, ['Looping over ' num2str(numel(datroi.trial)) ' trials'])
  
  for trl = 1:numel(datroi.trial);
    
    ft_progress(trl / numel(datroi.trial))
    x  = datroi.trial{trl}(r,:);
    
    %-----------------%
    %- 1. less than negthr
    npeaks = find(x < cfg.negthr);
    %-----------------%
    
    %-----------------%
    %- 2. distance between zero-crossing
    zcr = [0 diff(sign(x))]; % find all zero crossing
    
    %-------%
    %-negative peaks with zero crossing
    % np_zcr contains all the info for all the npeaks. Afterwards we discard
    % all the negative peaks that are part of the same slow wave (i.e. have
    % common zero-crossings
    np_zcr = nan(numel(npeaks), 5);
    %-------%
    
    for i = 1:numel(npeaks)
      
      %-------%
      %-previous zero-crossing
      zdwn = find( zcr(1:npeaks(i)) < 0, 1, 'last'); % preceding zero crossing
      zup = npeaks(i) + find( zcr(npeaks(i):end) > 0, 1, 'first'); % following zero crossing
      %-------%
      
      %-------%
      %-distance between crossing is between cfg.zcr (default .3 1s)
      if ~isempty(zup) && ~isempty(zdwn) && ...
          (zup - zdwn)/data.fsample >= cfg.zcr(1) && ...
          (zup - zdwn)/data.fsample <= cfg.zcr(2)
        
        np_zcr(i, :) = [i npeaks(i) x(npeaks(i)) zdwn zup];
      end
      %-------%
      
    end
    
    np_zcr(isnan(np_zcr(:,1)), :) = []; % delete empty rows at the end
    %-----------------%
    
    %-----------------%
    %-clean up duplicate peaks and call them slow waves
    [~, ~, unegw] = unique(np_zcr(:,4:5), 'rows'); % unique negative wave (check which peaks are in common to one slow wave)
    
    for sw = 1:numel( unique( unegw))
      swcnt = swcnt + 1;
      SW(swcnt).trl = trl;
      
      sw_zcr = np_zcr(unegw == sw, :);
      
      [neg_v, neg_i] = min(sw_zcr(:,3));
      
      SW(swcnt).negpeak_val = neg_v;
      SW(swcnt).negpeak_itrl  = sw_zcr(neg_i,2);
      SW(swcnt).negpeak_iabs  = datroi.sampleinfo(trl, 1) + SW(swcnt).negpeak_itrl;
      SW(swcnt).negpeak_time  = datroi.time{trl}(SW(swcnt).negpeak_itrl);
      
      SW(swcnt).begsw_itrl = sw_zcr(neg_i, 4);
      SW(swcnt).begsw_iabs = datroi.sampleinfo(trl, 1) + SW(swcnt).begsw_itrl;
      SW(swcnt).begsw_time = datroi.time{trl}(SW(swcnt).begsw_itrl);
      
      SW(swcnt).zcr_itrl = sw_zcr(neg_i, 5) -1;
      SW(swcnt).zcr_iabs = datroi.sampleinfo(trl, 1) + SW(swcnt).zcr_itrl -1;
      SW(swcnt).zcr_time = datroi.time{trl}(SW(swcnt).zcr_itrl);
      
    end
    %-----------------%
  end
  ft_progress('close')
  
end

if isempty(SW)
  return
end

%-----------------%
%-remove duplicates 
%-------%
%-pure duplicates
% (because of padding the same data is used in two consecutive trials, for example)
% check if the values of two negpeak in absolute samples are the same
dupl = [true diff([SW.negpeak_iabs]) ~= 0];
SW = SW(dupl);
%-------%

%-------%
%-similar, from different ROIs
[~, SWi] = sort([SW.negpeak_iabs]);
SW = SW(SWi);

% don't use only negpeak
dupl = [true diff([SW.negpeak_iabs]) > data.fsample/10];
SW = SW(dupl);
%-------%
%-----------------%

%---------------------------%

%---------------------------%
%- 3. peak-to-peak
%-----------------%
%-compute positive peak anyway
for sw = 1:numel(SW)
  
  %--------%
  %-get data
  begsw = SW(sw).begsw_itrl;
  endsw = SW(sw).zcr_itrl + cfg.postzcr * data.fsample;
  
  x  = mean(datroi.trial{SW(sw).trl},1);
  if endsw > size(x,2)
    endsw = size(x,2);
  end
  xsw = x(begsw:endsw);
  %--------%
  
  [pos_v, pos_i] = max(xsw);
  SW(sw).peak2peak = pos_v - SW(sw).negpeak_val;
  
  SW(sw).pospeak_val = pos_v;
  SW(sw).pospeak_itrl = SW(sw).begsw_itrl + pos_i - 1;
  SW(sw).pospeak_iabs = SW(sw).begsw_iabs + pos_i - 1;
  SW(sw).pospeak_time = datroi.time{SW(sw).trl}( SW(sw).pospeak_itrl );
end
%-----------------%

%-----------------%
if ~isempty(cfg.p2p)
    SWkeep = ([SW.pospeak_val] - [SW.negpeak_val]) >= cfg.p2p;
    SW = SW(SWkeep);
end
%-----------------%
%---------------------------%

%---------------------------%
%-EXTRA: traveling of slow waves
for i = 1:numel(SW)
  
  %-----------------%
  %-get the data
  datsel = SW(i).negpeak_itrl + cfg.trvl.wnw * data.fsample;

  endtrl = numel(data.time{SW(i).trl});
  if datsel(1) < 1; datsel(1) = 1; end
  if datsel(2) > endtrl;  datsel(2) = endtrl; end
  
  x = data.trial{SW(i).trl}(:, datsel(1):datsel(2));
  %-----------------%
  
  %-----------------%
  %-compute minimum
  [min_v, min_i] = min(x,[],2);
  
  %-------%
  %-keep channels above trvlnegthr
  swchan = min_v <= cfg.trvl.negthr;
  channel = find(swchan);
  label   = data.label(swchan);
  %-------%
  
  %-------%
  %-calculate delay
  ichan = min_i(swchan); % sample index of the min peak for all the channels with slow wave
  [delay_sample, delay_order] = sort(ichan);
  delay_sample = delay_sample - delay_sample(1); % reref to first electrode
  delay_second = delay_sample / data.fsample;
  
  %-sort according to delay
  channel = channel(delay_order);
  label   = label(delay_order);
  %-------%
  
  %-------%
  %-if the negative peak is too far between two electrodes, you only
  %consider the first part
  lastelec = find(diff(delay_second) > cfg.trvl.dist, 1);
  if ~isempty(lastelec)
    delay_second = delay_second(1:lastelec);
    channel = channel(1:lastelec);
    label = label(1:lastelec);
  end
  %-------%
  
  %-------%
  %-prepare structure
  SW(i).channel = channel;
  SW(i).label = label;
  SW(i).delay = delay_second;
  SW(i).maxdelay = max(delay_second);
  %-------%
  %-----------------%
  
end
%---------------------------%
