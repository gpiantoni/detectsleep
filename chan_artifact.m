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
