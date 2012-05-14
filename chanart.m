function [data, output] = chanart(cfg, data)
%CHANART reject and repair single channels
%
% CFG
%  .sleepepochs.chanart.auto(1).met: method to automatically find bad
%                                    channels ('var' 'range' 'diff')
%  .sleepepochs.chanart.auto(1).thr: in microvolts (10000, 3000, 1000)
%
%  .sens.file: file with EEG sensors. It can be sfp or mat
%  .sens.dist: distance between sensors to consider them neighbors (in the units of cfg.sens.file)
%
% IN
%   sleep data with all epochs
%
% OUT
%   sleep data with all epochs, but channels have been repaired

%--------------------------%
%-reject channels with difference methods
output = '';

for a = 1:numel(cfg.sleepepochs.chanart.auto)
  
  %------------------%
  %-automatic (above threshold in variance)
  switch cfg.sleepepochs.chanart.auto(a).met
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
  i_bad = find(allchan > cfg.sleepepochs.chanart.auto(a).thr);
  badchan{a} = data.label(i_bad);
  %-------%
  
  if ~isempty(badchan{a})
    %-------%
    %-output (sort bad channels depending on values of allchan)
    [~, s_bad] = sort(allchan(i_bad), 'descend');
    badname = '';
    for b = s_bad'
      badname = sprintf('%s %s (% 6d)', badname, badchan{a}{b}, allchan(i_bad(b)));
    end
    
    outtmp = sprintf('    %s (% 6d): min % 5.2f, median % 5.2f, mean % 5.2f, std % 5.2f, max % 5.2f\n    %g channels were bad: %s\n\n', ...
      cfg.sleepepochs.chanart.auto(a).met, cfg.sleepepochs.chanart.auto(a).thr, min(allchan), median(allchan), mean(allchan), std(allchan), max(allchan), ...
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
%-repair channels
%-------%
%-create neighbors from file
sens = ft_read_sens(cfg.sens.file);
sens.label = upper(sens.label);

cfg1 = [];
cfg1.elec = sens;
cfg1.method = 'distance';
cfg1.neighbourdist = cfg.sens.dist;
[~, neigh] = evalc('ft_prepare_neighbours(cfg1);');
%-------%

allbad = cat(1, badchan{:});
cfg1 = [];
cfg1.badchannel = allbad;
cfg1.neighbours = neigh;
data.elec = sens;
[~, data] = evalc('ft_channelrepair(cfg1, data);');
data = rmfield(data, 'elec');
%-----------------%
