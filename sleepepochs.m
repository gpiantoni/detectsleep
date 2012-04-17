function sleepanalysis(cfg, subj)
%SLEEPANALYSIS read single subject sleep data and apply analysis tools
% cfg.data, cfg.mod, cfg.proj, cfg.rec,
% cfg.stage, cfg.scorer

%---------------------------%
%-start log
output = sprintf('(p%02.f) %s started at %s on %s\n', ...
  subj,  mversion, datestr(now, 'HH:MM:SS'), datestr(now, 'dd-mmm-yy'));
tic_t = tic;
%---------------------------%

%---------------------------%
%-dir and files
cond = 'sleep';
ddir = sprintf('%s%04.f/%s/%s/', cfg.data, subj, cfg.mod, cond); % data
dfile = sprintf('%s_%s_%04.f_%s_%s', cfg.proj, cfg.rec, subj, cfg.mod, cond);
%---------------------------%

%-----------------------------------------------%
%-loop over sleep stage
for s = cfg.stage
  
  %---------------------------%
  %-read one epoch at the time
  D = spm_eeg_load([ddir dfile]);
  
  score = D.CRC.score{1, cfg.scorer};
  epch = find(ismember(score, s));
  
  %-------%
  %-TODO: feedback
  % using scorer, number of epoch in this stage
  %-------%
  
  %-------%
  %-find artifacts
  if ~isempty(D.CRC.score{5,cfg.scorer})
    artbeg = round(D.CRC.score{5,cfg.scorer}(:,1) * fsample(D)); % from time into samples
    artend = round(D.CRC.score{5,cfg.scorer}(:,2) * fsample(D)); % from time into samples
  else
    output = sprintf('%sWARNING: is the artifact rejection of scorer 3 empty???\n', output);
    artbeg = round(D.CRC.score{5,1}(:,1) * fsample(D)); % from time into samples
    artend = round(D.CRC.score{5,1}(:,2) * fsample(D)); % from time into samples
  end
  art = [artbeg artend];
  %-------%
  %---------------------------%
  
  %-------------------------------------%
  %-loop over epch
  swall = [];
  spall = [];
  
  %-----------------%
  %-try to estimate the size of the freqall variable
  %TODO: preallocation (especially for freq analysis)
  %-----------------%
  
  cnt = 0; % epoch count
  ngood = 0; % it's used to calculate the mean just before saving
  nbad = 0; % bad electrodes
  
  for e = epch
    
    %-------%
    %-progress
    cnt = cnt + 1;
    if strcmp(cfg.detsw.feedback, 'yes')
      fprintf('% 4.f/% 4.f\n', cnt, numel(epch))
    end
    %-------%
    
    %---------------------------%
    %-convert data and preprocessing
    %-----------------%
    %-convert
    cfg1 = [];
    cfg1.epoch = e;
    cfg1.pad = cfg.pad;
    data = sleep2ft(cfg1, [ddir dfile]);
    data.label = layout.label(4:end-2);
    %-----------------%
    
    %-----------------%
    %-reject if there is artifact
    if ~strcmp(cfg.reject, 'no')
      try
        cfg2 = [];
        cfg2.artfctdef.manual.artifact = art;
        cfg2.artfctdef.reject = cfg.reject;
        cfg2.artfctdef.minaccepttim = cfg.freqsw.length;
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
    %---------------------------%
    
    %---------------------------%
    %-reject bad sleep
    %-----------------%
    %-clean bad channels
    [data outtmp] = artchan(cfg, data);
    if ~isempty(outtmp)
      % output = sprintf('%se% 4.f%s', output, e, outtmp);
      nbad = nbad + outtmp;
    end
    %-----------------%
    %---------------------------%
    
    %-----------------%
    %-reref
    cfg2 = cfg.reref;
    [~, data] = evalc('ft_preprocessing(cfg2, data)');
    %-----------------%
    
    %-----------------%
    %-detect slow waves
    if any(strcmp(cfg.rundet, 'sw'))
      cfg3 = cfg.detsw.sw;
      [sw] = sw_detect(cfg3, data);
      
      if ~isempty(sw)
        [sw.trl] = deal(e);
        swall = [swall sw];
      end
    end
    %-----------------%
    
    %-----------------%
    %-detect spindles
    if any(strcmp(cfg.rundet, 'sp'))
      cfg4 = cfg.detsw.sp;
      cfg4.roi = cfg.detsw.sw.roi;
      [sp] = sp_detect(cfg4, data);
      
      if ~isempty(sp)
        [sp.trl] = deal(e);
        spall = [spall sp];
      end
    end
    %-----------------%
    
    %-----------------%
    %-freq analysis
    %-another smart way to go about this is to delete parts first and then
    %use redefinetrial
    if any(strcmp(cfg.rundet, 'freq'))
      
      %-------%
      %-split in 2-s time windows
      cfg5 = [];
      cfg5.length = cfg.freqsw.length;
      cfg5.overlap = .5;
      data = ft_redefinetrial(cfg5, data);
      %-------%
      
      %-------%
      %-goodtrl, exclude the last one
      %in this way, it uses the 30s epoch (the first epoch is -1:1, so you have
      %one second of the previous epoch, even if the previous epoch is not a
      %selected sleep state. There are so few of these epochs and it's not a
      %problem)
      if cfg.pad == 1
        goodtrl = 1: numel(data.trial)-1;
      else
        error('goodtrl is undefined if cfg.pad ~= 1')
      end
      %-------%
      
      %       %-------%
      %       %-find good trl
      %       if strcmp(cfg.reject, 'partial')
      %
      %         goodlog = true(size(goodtrl));
      %         for t = goodtrl
      %
      %           if hasart(data.time{t}(1), data.time{t}(end), artbeg, artend)
      %             goodlog(t) = false;
      %
      %           end
      %         end
      %
      %         if ~any(goodlog)
      %           outtmp = sprint('the whole epoch % 4.f was rejected (THIS SHOULD NOT HAPPEN!!!)\n', e);
      %           output = [output outtmp];
      %           continue
      %
      %         elseif ~all(goodlog) % not all are good
      %           outtmp = sprintf('% 3.f segments of epoch % 4.f were rejected\n', numel(find(goodlog == false)), e);
      %           output = [output outtmp];
      %           goodtrl = goodtrl(goodlog);
      %         end
      %
      %       end
      %       %-------%
      
      for r = 1:numel(cfg.freqsw.roi)
        
        %-------%
        cfg6 = [];
        cfg6.method = 'mtmfft';
        cfg6.foilim = cfg.freqsw.foilim;
        cfg6.taper = 'hanning';
        cfg6.feedback = 'none';
        cfg6.channel = cfg.freqsw.roi(r).chan;
        cfg6.trials = goodtrl;
        freq = ft_freqanalysis(cfg6, data);
        %-------%
        
        freqall{r}(cnt, :) = squeeze(mean(freq.powspctrm,1)) * numel(goodtrl);
        ngood = ngood + numel(goodtrl);
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
  if any(strcmp(cfg.rundet, 'sw'))
    %-------%
    %-pure duplicates
    % (because of padding the same data is used in two consecutive trials, for example)
    % check if the values of two negpeak in absolute samples are the same
    dupl = [true diff([swall.negpeak_iabs]) ~= 0];
    swall = swall(dupl);
    sw = swall;
    %-------%
    
    %-------%
    %-remove sw during artifacts
    if strcmp(cfg.reject, 'partial')
      artsw = false(numel(sw),1);
      for i = 1:numel(artbeg)
        artsw = artsw | artbeg(i) < [sw.negpeak_iabs]' & [sw.negpeak_iabs]' < artend(i);
      end
      
      if numel(find(artsw))
        artswstr = sprintf('%10.f ', [sw(artsw).trl]);
        outtmp = sprintf('ERROR: % 3.f slow waves in bad segments out of % 5.f, at trials: %s \n', ...
          numel(find(artsw)), numel(sw), artswstr);
        output = [output outtmp];
      end
      
      sw(artsw) = [];
    end
    %-------%
    
    %-------%
    %-save
    SWfile = sprintf('%sdetSW_%1.f_%03.f', cfg.detd, s, subj);
    save(SWfile, 'sw')
    %-------%
  end
  %-----------------%
  
  %-----------------%
  %-spindles
  if any(strcmp(cfg.rundet, 'sp'))
    %-------%
    %-pure duplicates
    % (because of padding the same data is used in two consecutive trials, for example)
    % check if the values of two negpeak in absolute samples are the same
    dupl = [true diff([spall.maxsp_iabs]) ~= 0];
    spall = spall(dupl);
    sp = spall;
    %-------%
    
    %-------%
    %-remove sw during artifacts
    if strcmp(cfg.reject, 'partial')
      artsp = false(numel(sp),1);
      for i = 1:numel(artbeg)
        artsp = artsp | artbeg(i) < [sp.maxsp_iabs]' & [sp.maxsp_iabs]' < artend(i);
      end
      
      if numel(find(artsp))
        artspstr = sprintf('%10.f ', [sp(artsp).trl]);
        outtmp = sprintf('ERROR: % 3.f spindles in bad segments out of % 5.f, at samples: %s \n', ...
          numel(find(artsp)), numel(sp), artspstr);
        output = [output outtmp];
      end
      
      sp(artsp) = [];
    end
    %-------%
    
    %-------%
    %-save
    SPfile = sprintf('%sdetSP_%1.f_%03.f', cfg.detd, s, subj);
    save(SPfile, 'sp')
    %-------%
  end
  %-----------------%
  
  %-----------------%
  %-average over epoch
  if any(strcmp(cfg.rundet, 'freq'))
    
    freqfile = sprintf('%sfreq_%1.f_%03.f', cfg.detd, s, subj);
    for r = 1:numel(cfg.freqsw.roi)
      freqname = ['freq_' cfg.freqsw.roi(r).name];
      eval([freqname ' = mean(freqall{r}, 1) / ngood;'])
      
      if ~exist([freqfile '.mat'], 'file')
        save(freqfile, freqname, 'ngood')
      else
        save(freqfile, freqname, '-append')
      end
      
    end
  end
  %-----------------%
  
  output = sprintf('%saverage number of bad channels: % 5.f\n', output, nbad / numel(epch));
  %---------------------------%
  
  clear swall sw spall sp freq freqall ngood
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

function [arttype] = hasart(trlbeg, trlend, artbeg, artend)
%HASART check if there are artifacts in the trial

arttype = 0;

if any(trlbeg < artbeg & trlend > artbeg)
  arttype = 1;
end

if any(trlbeg < artend &trlend > artend)
  arttype = 2;
end

if any(trlbeg > artbeg & trlend < artend)
  arttype = 3;
end

function [data, output] = artchan(cfg, data)

%--------------------------%
%-reject channels with difference methods
output = [];

for a = 1:numel(cfg.badchan.auto)
  
  %------------------%
  %-automatic (above threshold in variance)
  switch cfg.badchan.auto(a).met
    case 'var'
      %-------%
      %-compute var
      allchan = std([data.trial{:}], [], 2).^2;
      %-------%
      
    case 'range'
      %-------%
      %-compute range
      alldat = [data.trial{:}];
      allchan = range(alldat,2);
      %-------%
      
    case 'diff'
      %-------%
      %-compute range
      alldat = [data.trial{:}];
      allchan = max(abs(diff(alldat')))';
      %-------%
      
  end
  %------------------%
  
  %------------------%
  %-find badchan and output
  %-------%
  %-define bad channels
  i_bad = feval(eval(cfg.badchan.auto(a).fun), allchan);
  badchan{a} = data.label(i_bad);
  %-------%
  
  if ~isempty(badchan{a})
    %-------%
    %-output (sort bad channels depending on values of allchan)
    [~, s_bad] = sort(allchan(i_bad), 'descend');
    badname = '';
    for b = s_bad'
      badname = sprintf('%s %s (%6.f)', badname, badchan{a}{b}, allchan(i_bad(b)));
    end
    
    outtmp = sprintf('    %s (%s): min % 5.2f, median % 5.2f, mean % 5.2f, std % 5.2f, max % 5.2f\n    %g channels were bad: %s\n\n', ...
      cfg.badchan.auto(a).met, cfg.badchan.auto(a).fun, min(allchan), median(allchan), mean(allchan), std(allchan), max(allchan), ...
      numel(badchan{a}), badname);
    output = [output outtmp];
  end
  %-------%
  %------------------%
  
end
%--------------------------%

%-----------------%
%-do not repair channels 
if isempty(output)
  return
end
%-----------------%

%-----------------%
%-all bad and check if the reference is bad
allbad = cat(1, badchan{:});
output = numel(unique(allbad));

% if ~any(strcmp('all', cfg.reref.refchannel))
%   
%   badref = intersect(cfg.reref.refchannel, allbad);
%   
%   if ~isempty(badref)
%     badrefstr = sprintf(' %s,', badref{:});
%     outtmp = sprintf('  WARNING: Reference channel (%s) is bad\n', badrefstr);
%     output = [output outtmp];
%   end
%   
% end
%-----------------%

%-----------------%
%-repair channels
load(cfg.elecfile, 'elec', 'nbor')
cfg1 = [];
cfg1.badchannel = allbad;
cfg1.neighbours = nbor;
data.elec = elec;
[data] = ft_channelrepair(cfg1, data);
data = rmfield(data, 'elec');
%-----------------%
