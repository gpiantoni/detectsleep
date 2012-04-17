function data = sleep2ft(cfg, D)
%SLEEP2FT read the data in FASST format into FieldTrip
%
% data = sleep2ft(cfg, D)
%
% cfg
%   .stage = [3 4]
%   .scorer = 1 (string or index of D.CRC.score)
%  OR
%   .epoch = index of the epoch(s) to select
%
%  (optional)
%   .pad = amount of padding (s)
%   .chantype = 'EEG' 
%   (or .channel  = 'e016' or index of channels or cell with labels)
%   .cat = concatenate into one large file ('yes' or 'no')
%
% D: can be string or MEEG/D
%
% data: in fieldtrip format
%   .trialinfo contains epoch index and epoch scoring when more trials or
%   only scoring when using concatenation

%---------------------------%
%-input check
%-----------------%
%-cfg
if ~isfield(cfg, 'pad'); cfg.pad = 0; end
if ~isfield(cfg, 'cat'); cfg.cat = 'no'; end
if cfg.pad ~= 0 && strcmp(cfg.cat, 'yes')
  error('it doesn''t make sense to add padding and concatenate the epochs')
end
%-----------------%

%-----------------%
%-D
if ischar(D)
  D = spm_eeg_load(D);
elseif isstruct(D)
  D = meeg(D);
end

if ~isfield(D, 'CRC')
  error(['Sleep file does not have CRC field'])
end
%-----------------%
%---------------------------%

%---------------------------%
%-find epochs
if isfield(cfg, 'epoch')
  epch = cfg.epoch;
  
  %-----------------%
  %-try to use indicative scorer, otherwise NaN
  if isfield(D.CRC, 'score')
    score = D.CRC.score{1,1};
  else
    error('no scoring in the files')
    % ... create NaN: not strictly necessary, only for trialinfo
  end
  %-----------------%
else
  %-----------------%
  %-use scorer and stage from cfg
  %-------%
  %-find scorer as index
  if ischar(cfg.scorer)
    scorer = find(strcmp(D.CRC.score(2,:), cfg.scorer));
    
    if numel(scorer) ~= 1
      error(['could not find ' cfg.scorer ' in D.CRC.score'])
    end
    
  else
    scorer = cfg.scorer;
  end
  %-------%
  
  score = D.CRC.score{1, scorer};
  epch = find(ismember(score, cfg.stage));
  %-----------------%
  
end
%---------------------------%

%---------------------------%
%-create FT structure
data = [];
data.fsample = fsample(D);

%-----------------%
%-labels
data.label = chanlabels(D)';
if isfield(cfg, 'chantype') && ischar(cfg.chantype)
  chan = find(strcmp(chantype(D), cfg.channel));
  data.label = data.label(chan);
  
elseif isfield(cfg, 'channel') && ischar(cfg.channel)
  chan = find(strcmp(chanlabels(D), cfg.channel));
  data.label = data.label(chan);

elseif isfield(cfg, 'channel') && iscellstr(cfg.channel)
  [~, chan] = intersect(data.label, cfg.channel);
  data.label = data.label(chan);

elseif isfield(cfg, 'channel')
  chan = cfg.channel;
  data.label = data.label(chan);

else
  chan = find(true(size(data.label)));
end
%-----------------%

%-----------------%
%-epochs
for i = 1:numel(epch)
  e = epch(i);
  
  %------%
  %-time
  begtime = (e-1) * D.CRC.score{3,1} - cfg.pad + 1/fsample(D);
  endtime = e * D.CRC.score{3,1} + cfg.pad;
  time{i} = begtime:1/fsample(D):endtime;
  %------%
  
  %------%
  %-trial
  begsam = begtime * fsample(D);
  endsam = endtime * fsample(D);
  trial{i} = D(chan, begsam:endsam,1);
  %------%
  
  %------%
  %-extra info into
  sampleinfo(i, 1:2) = round([begsam endsam]); % avoid weird rounding errors
  trialinfo(i, 1:2) = [e score(e)];
  %------%
  
end
%-----------------%

%-----------------%
%-concatenate or not
if strcmp(cfg.cat, 'yes')
  
  data.time{1} = [time{:}];
  data.trial{1} = [trial{:}];
  data.trialinfo = trialinfo(:,2)'; % keep only scoring value (transpose bc only one trial)
  
else
  
  data.time = time;
  data.trial = trial;
  data.trialinfo = trialinfo;
  data.sampleinfo = sampleinfo;
  
end
%-----------------%
%---------------------------%