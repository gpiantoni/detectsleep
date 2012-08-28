function data = fasst2fieldtrip(cfg, D)
%FASST2FIELDTRIP read the data in FASST format into FieldTrip
%
% data = sleep2ft(cfg, D)
%
% cfg
%   .stage = [3 4]
%   .rater = 1 (string or index of D.CRC.score)
%  OR
%   .epoch = index of the epoch(s) to select
%
%  (optional)
%   .preproc: any option that you want to pass to ft_preprocessing
%
% D: can be string to FASST file or MEEG/D
%
% data: in fieldtrip format
%   .trialinfo contains epoch index and epoch scoring when more trials or
%   only scoring when using concatenation

%---------------------------%
%-input check
%-----------------%
%-cfg
if ~isfield(cfg, 'preproc'); cfg.preproc = []; end
%-----------------%

%-----------------%
%-D
if ischar(D)
  load(D, 'D');
end

D = struct(D); % don't use OO of SPM8

if ~isfield(D.other, 'CRC')
  error(['Sleep file does not have CRC field'])
else
  score = D.other.CRC.score;
end
%-----------------%
%---------------------------%

%---------------------------%
%-find epochs
if isfield(cfg, 'epoch')
  epch = cfg.epoch;
  
else
  
  %-----------------%
  %-use rater and stage from cfg
  %-------%
  %-find rater as index
  if ischar(cfg.rater)
    rater = find(strcmp(score(2,:), cfg.rater));
    
    if numel(rater) ~= 1
      error(['could not find ' cfg.rater ' in D.CRC.score'])
    end
    
  else
    rater = cfg.rater;
  end
  %-------%
  
  score = score(:, rater);
  epch = find(ismember(score{1}, cfg.stage));
  %-----------------%
  
end
%---------------------------%

%---------------------------%
%-read data
%-----------------%
% epoch has to  be vertical
if size(epch, 1) == 1
  epch = epch';
end
%-----------------%

%-----------------%
%-time info
wndw = score{3} * D.Fsample;
beginsleep = score{4}(1) * D.Fsample;
begsample = (epch - 1) * wndw + beginsleep;
endsample = begsample + wndw - 1;

trl = [round([begsample endsample begsample]) epch score{1}(epch)'];
%-----------------%

%-----------------%
tmpcfg = cfg.preproc;
tmpcfg.dataset = [D.path D.fname];
tmpcfg.trl = trl;
data = ft_preprocessing(tmpcfg);
%---------------------------%