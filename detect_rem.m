function [rem] = detect_rem(cfg, data)
%DETECT_REM detect REM using channel electrodes
% Simple preprocessing. I suggest at the moment to do 
%   1_ low pass filtering
%   2_ take first derivative
%   3_ detect REM above a threshold
%
% Use as:
%    [rem] = detect_rem(cfg, data)
%
% cfg
%  two channel groups for bipolar montage
%  .eog(1).name: 'eog'
%  .eog(1).eog1: cell with labels (better)
%  .eog(1).eog2: cell with labels (better)
%
%  .preproc: cfg to pass to ft_preprocessing
%  now, it includes:
%  .preproc.lpfilter = 'yes';
%  .preproc.lpfreq = 15; % Hz
%  .preproc.derivative = 'yes';
%  
%  .thr: threshold to include REM
%
% data
%    data in fieldtrip format, see also sleep2ft
%
% REM
%    struct for each detected rem
%  .trl: trial it belongs to
%
%  .label: which label was used for this detected REM
%
%  .begin_itrl: first point above threshold in samples, from beginning of the trial
%  .begin_iabs: first point above threshold in samples, based on sampleinfo
%  .begin_time: first point above threshold in seconds, based on data.time
%
%  .end_itrl: first point above threshold in samples, from beginning of the trial
%  .end_iabs: first point above threshold in samples, based on sampleinfo
%  .end_time: first point above threshold in seconds, based on data.time
%

%---------------------------%
%-prepare input
%-----------------%
%-check cfg
%-------%
%-defaults
if ~isfield(cfg, 'preproc'); cfg.preproc = []; end
if ~isfield(cfg.preproc, 'lpfilter'); cfg.preproc.lpfilter = 'yes'; end
if ~isfield(cfg.preproc, 'lpfreq'); cfg.preproc.lpfreq = 20; end
if ~isfield(cfg.preproc, 'derivative'); cfg.preproc.derivative = 'yes'; end
if ~isfield(cfg, 'thr'); cfg.thr = 15; end
%-------%
%-----------------%
%---------------------------%

%---------------------------%
%-prepare bipolar montage
mont = [];
mont.labelorg = data.label;
mont.labelnew = {cfg.eog.name};
mont.tra = zeros(numel(mont.labelnew), numel(mont.labelorg));

for i = 1:numel(cfg.eog)
  
  [~, eog1] = intersect(mont.labelorg, cfg.eog(i).eog1);
  [~, eog2] = intersect(mont.labelorg, cfg.eog(i).eog2);
  
  mont.tra(i, eog1) = 1;
  mont.tra(i, eog2) = -1;
  
end

data = ft_apply_montage(data, mont);
%---------------------------%

%---------------------------%
%-filtering
cfg1 = cfg.preproc;
data = ft_preprocessing(cfg1, data);
%---------------------------%

%---------------------------%
%-REM is above threshold
rem = [];
cnt = 0;

for t = 1:numel(data.trial)
  
  for i = 1:numel(data.label)
    
    x = data.trial{t}(i,:);
    
    %-----------------%
    %-above threshold
    abrem = find(abs(x) >= cfg.thr);
    if numel(abrem) ==0
      continue
    end
    %-----------------%
    
    %-----------------%
    %-collect sample info
    begrem = [1 find([1 diff(abrem)] ~= 1)];
    endrem = [find([1 diff(abrem)] ~= 1)-1 numel(abrem)];
    %-----------------%
    
    for r = 1:numel(begrem)
      cnt = cnt + 1;
      
      %-----------------%
      %-prepare output variable
      rem(cnt).trl = t;
      rem(cnt).begin_itrl = abrem(begrem(r));
      rem(cnt).begin_iabs = abrem(begrem(r)) + data.sampleinfo(t, 1) - 1;
      rem(cnt).begin_time = data.time{t}(abrem(begrem(r)));
      
      rem(cnt).end_itrl = abrem(endrem(r));
      rem(cnt).end_iabs = abrem(endrem(r)) + data.sampleinfo(t, 1) - 1;
      rem(cnt).end_time = data.time{t}(abrem(endrem(r)));
      %-----------------%
      
    end
  end
  
end
%---------------------------%
