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
  %-use scorer and stage from cfg
  %-------%
  %-find scorer as index
  if ischar(cfg.scorer)
    scorer = find(strcmp(score(2,:), cfg.scorer));
    
    if numel(scorer) ~= 1
      error(['could not find ' cfg.scorer ' in D.CRC.score'])
    end
    
  else
    scorer = cfg.scorer;
  end
  %-------%
  
  score = score(1, scorer);
  epch = find(ismember(score, cfg.stage));
  %-----------------%
  
end
%---------------------------%

%---------------------------%
%-read data
% epoch has to  be vertical
%-----------------%
%-time info
wndw = score{3,1} * D.Fsample;
beginsleep = score{4,1}(1) * D.Fsample;
begsample = (epoch - 1) * wndw + beginsleep;
endsample = begsample + wndw - 1;

trl = [round([begsample endsample]) zeros(numel(epoch),1) epoch score{1}(epoch)'];
%-----------------%

%-----------------%
tmpcfg = cfg.preproc;
tmpcfg.dataset = [D.path D.fname];
tmpcfg.trl = trl;
data = ft_preprocessing(tmpcfg);
%---------------------------%