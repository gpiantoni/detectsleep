function sleepepochs(cfg, subj)
%SLEEPEPOCHS analyze sleep data epoch by epoch
%
% CFG
%  .rec: name of the recording
%  .data: name of projects/PROJNAME/subjects/
%  .mod: recording modality (likely 'eeg')
%
%  .stage: which stages should be analyzed (REM = 5)
%  .scorer: index of the scorer, for the stages
%
%  .sleepepochs.feedback: gives feedback or no of each epoch (logical, default: false)
%  .sleepepochs.pad: amount of padding for each epoch (default: 1s)
%  .sleepepochs.visrej: use visual rejection from sleep scoring (logical)
%  .sleepepochs.chanrej: reject channels outside the limit (TODO)
%
%  .sleepepochs.reref: if not empty, re-referencing
%  .sleepepochs.reref.refchannel: channels used for re-referencing ({'E94' 'E190'} or 'all')
%  .sleepepochs.reref.implicit: implicit channel ('E257')
%
%  .sleepepochs.detdir: directory where you want to save the slow waves
%  .sleepepochs.detsw: if not empty, slow wave detection, you should use the
%                      configuration options based on DETECT_SLOWWAVE
%  .sleepepochs.detsp: if not empty, spindle detection, you should use the
%                      configuration options based on DETECT_SPINDLE
%
% INPUT
%  - single-subject sleep data, in the format gosd_svui_XXXX_eeg_sleep.mat
%    but in the future, the projname (GOSD) will probably go.
%
% OUTPUT
%  - depeding on the detection routine, you'll get in .sleepepochs.detdir:
%    slowwave_stageX_XXXX.mat with slow waves for each subject
%    spindle_stageX_XXXX.mat with spindles for each subject
%
% Part of DETECTSLEEP
% see also SLEEPEPOCHS, SLEEP2FT, DETECT_ARTIFACT
%          DETECT_SPINDLE, DETECT_SLOWWAVE, SPINDLE_PEAK, SLEEP_FREQ

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s started at %s on %s\n', ...
  subj, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-default parameters
if ~isfield(cfg, 'sleepepochs'); cfg.sleepepochs = []; end
if ~isfield(cfg.sleepepochs, 'feedback'); cfg.sleepepochs.feedback = false; end
if ~isfield(cfg.sleepepochs, 'pad'); cfg.sleepepochs.pad = 1; end
if ~isfield(cfg.sleepepochs, 'visrej'); cfg.sleepepochs.visrej = true; end
if ~isfield(cfg.sleepepochs, 'chanrej'); cfg.sleepepochs.chanrej = false; end
%---------------------------%

%---------------------------%
%-dir and files
cond = 'sleep';
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cond); % data
dfile = sprintf('%s_%s_%04.f_%s_%s', 'gosd', cfg.rec, subj, cfg.mod, cond); % TODO: i think that the projname should go, but we need to check how the files were created

load(cfg.sens.layout, 'layout') % rename labels for consistency
%---------------------------%

%-----------------------------------------------%
%-loop over sleep stage
for s = cfg.stage
  
  %---------------------------%
  %-read one epoch at the time
  D = spm_eeg_load([ddir dfile]);
  
  score = D.CRC.score{1, cfg.scorer};
  epch = find(ismember(score, s));
  
  %-----------------%
  %-feedback
  outtmp = sprintf('scorer ''%s'', stage %d, number of epochs % 4d\n', ...
    D.CRC.score{2,cfg.scorer}, s, numel(epch));
  output = [output outtmp];
  %-----------------%
  
  %-----------------%
  %-read visually detected artifacts
  if ~isempty(D.CRC.score{5,cfg.scorer})
    artbeg = round(D.CRC.score{5,cfg.scorer}(:,1) * fsample(D)); % from time into samples
    artend = round(D.CRC.score{5,cfg.scorer}(:,2) * fsample(D)); % from time into samples
  else
    output = sprintf('%sWARNING: is the artifact rejection of scorer 3 empty???\n', output);
    artbeg = round(D.CRC.score{5,1}(:,1) * fsample(D)); % from time into samples
    artend = round(D.CRC.score{5,1}(:,2) * fsample(D)); % from time into samples
  end
  art = [artbeg artend];
  %-----------------%
  %---------------------------%
  
  %-------------------------------------%
  %-loop over epch
  %TODO: preallocation (especially for freq analysis)
  swall = [];
  
  cnt = 0; % epoch count
  ngood = 0; % it's used to calculate the mean just before saving
  nbad = 0; % bad electrodes
  
  for e = epch
    
    %-----------------%
    %-progress
    cnt = cnt + 1;
    if cfg.sleepepochs.feedback
      fprintf('stage %d/%d (% 4d/% 4d)\n', ...
        find(cfg.stage==s), numel(cfg.stage), cnt, numel(epch))
    end
    %-----------------%
    
    %---------------------------%
    %-convert data and preprocessing
    cfg1 = [];
    cfg1.epoch = e;
    cfg1.pad = cfg.sleepepochs.pad;
    data = sleep2ft(cfg1, [ddir dfile]);
    
    data.label = layout.label(4:end-2);
    %---------------------------%
    
    %---------------------------%
    %-reject artifacts
    %-----------------%
    %-visual rejection during visual scoring
    if cfg.sleepepochs.visrej
      try
        cfg2 = [];
        cfg2.artfctdef.manual.artifact = art;
        cfg2.artfctdef.reject = 'partial';
        cfg2.artfctdef.minaccepttim = 2;
        data = ft_rejectartifact(cfg2, data);
      catch
        output = sprintf('%sComplete rejection of epoch % 3.f\n', output, e);
        continue
      end
    end
    %-----------------%
    
    %-----------------%
    %-remove nan
    for i = 1:numel(data.trial)
      data.trial{i}(isnan(data.trial{i})) = 0;
    end
    %-----------------%
    
    %-----------------%
    %-clean bad channels
    if cfg.sleepepochs.chanrej
      [data outtmp] = chanart(cfg, data);
      output = [output outtmp];
    end
    %-----------------%
    %---------------------------%
    
    %---------------------------%
    %-reref
    if isfield(cfg.sleepepochs, 'reref') && ~isempty(cfg.sleepepochs.reref)
      cfg3 = cfg.sleepepochs.reref;
      [~, data] = evalc('ft_preprocessing(cfg3, data)');
    end
    %---------------------------%
    
    %---------------------------%
    %-detection
    %-----------------%
    %-detect slow waves
    if isfield(cfg.sleepepochs, 'detsw') && ~isempty(cfg.sleepepochs.detsw)
      [sw] = detect_slowwave(cfg.sleepepochs.detsw, data);
      
      if ~isempty(sw)
        [sw.trl] = deal(e);
        swall = [swall sw];
      end
      
    end
    %-----------------%
    
    %-----------------%
    %-detect spindles
    if isfield(cfg.sleepepochs, 'detsp') && ~isempty(cfg.sleepepochs.detsp)
      [sp] = detect_spindle(cfg.sleepepochs.detsp, data);
      
      if ~isempty(sp)
        [sp.trl] = deal(e);
        spall = [spall sp];
      end
    end
    %-----------------%
    %---------------------------%
    
  end
  %-------------------------------------%
  
  %---------------------------%
  %-save file
  %-----------------%
  %-slow waves
  if isfield(cfg.sleepepochs, 'detsw') && ~isempty(cfg.sleepepochs.detsw)
    
    %-------%
    %-pure duplicates
    % (because of padding the same data is used in two consecutive trials, for example)
    % check if the values of two negpeak in absolute samples are the same
    dupl = [true diff([swall.negpeak_iabs]) ~= 0];
    swall = swall(dupl);
    %-------%
    
    %-------%
    %-save
    slowwave = swall;
    swfile = sprintf('%sslowwave_stage%1d_%04d', cfg.sleepepochs.detsw.dir, s, subj);
    save(swfile, 'slowwave')
    %-------%
    
  end
  %-----------------%
  
  %-----------------%
  %-spindles
  if isfield(cfg.sleepepochs, 'detsp') && ~isempty(cfg.sleepepochs.detsp)
    
    %-------%
    %-pure duplicates
    % (because of padding the same data is used in two consecutive trials, for example)
    % check if the values of two negpeak in absolute samples are the same
    dupl = [true diff([spall.maxsp_iabs]) ~= 0];
    spall = spall(dupl);
    %-------%
    
    %-------%
    %-save
    spindle = spall;
    spfile = sprintf('%sspindle_stage%1d_%04d', cfg.sleepepochs.detsw.dir, s, subj);
    save(spfile, 'spindle')
    %-------%
  end
  %-----------------%
  
  output = sprintf('%saverage number of bad channels: % 5.f\n', output, nbad / numel(epch));
  %---------------------------%
  
  clear swall sw slowwave spall sp spindle
end

%---------------------------%
%-end log
toc_t = toc(tic_t);
outtmp = sprintf('(p%02.f) %s ended at %s on %s after %s\n\n', ...
  subj, mfilename, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'), ...
  datestr( datenum(0, 0, 0, 0, 0, toc_t), 'HH:MM:SS'));
output = [output outtmp];

%-----------------%
fprintf(output)
fid = fopen([cfg.log '.txt'], 'a');
fwrite(fid, output);
fclose(fid);
%-----------------%
%---------------------------%
