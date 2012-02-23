function [peak] = find_spindles_peakfreq(cfg, data)
%FIND_SPINDLES_PEAKFREQ find max power in spindle frequency band
% Use as:
%    [peak] = find_spindles_peakfreq(cfg, data)
% where
%  cfg has optional fields
%  .output = 'average' (default) or 'epochs' (one peak for the whole recording or one for each epoch)
%  .foilim = frequency to look for spindles (default [8 18])
%  .length = length of smaller windows (1/cfg.length is the frequency resolution, default 5)
%  .plot   = true or false, only for output average, it plots the power spectrum (default: false)
%  .order  = order of the model used to remove 1/f signal (default 2)

% 11/10/26 gp: created

%-------------------------------------%
%-check input
if ~isfield(cfg, 'length'); cfg.length = 5; end
if ~isfield(cfg, 'foilim'); cfg.foilim = [8 18]; end
if ~isfield(cfg, 'order'); cfg.order = 1; end
if ~isfield(cfg, 'output'); cfg.output = 'average'; end
if ~isfield(cfg, 'plot'); cfg.plot = false; end

%-----------------%
%-prepare average/epochs
feedback = 'none';
epochs = 1:numel(data.trial);
if strcmp(cfg.output, 'average')
  epochs = epochs';
  feedback = 'etf';
end

if strcmp(cfg.output, 'epochs')
  cfg.plot = false;
end
%-----------------%
%-------------------------------------%

%---------------------------------------------------------%
%-loop over epochs
peak = NaN(size(epochs,2), numel(data.label));

for e = epochs
  
  if strcmp(cfg.output, 'epochs')
    disp(e/numel(epochs)*100)
  end
  
  %-------------------------------------%
  cfg1 = [];
  cfg1.trials = e;
  cfg1.length = cfg.length;
  cfg1.overlap = .5;
  [~, shortdata] = evalc('ft_redefinetrial(cfg1, data);');
  %-------------------------------------%
  
  %-------------------------------------%
  cfg1 = [];
  cfg1.method = 'mtmfft';
  cfg1.taper = 'hanning';
  cfg1.foilim = cfg.foilim;
  cfg1.feedback = feedback;
  if strcmp(feedback, 'etf')
    freq = ft_freqanalysis(cfg1, shortdata);
  else % hide feedback
    [~, freq] = evalc('ft_freqanalysis(cfg1, shortdata);');
  end
  %-------------------------------------%
  
  %-------------------------------------%
  %-loop over channels
  for c = 1:numel(data.label)
    b = log( freq.powspctrm(c, :));
    
    %-----------------%
    %-removing linear trends with GLM (especially good for 1/f noise)
    basis    = 1:numel(b);
    x        = zeros(cfg.order+1, numel(b));
    for i=0:cfg.order
      x(i+1,:) = basis.^(i);
    end
    b = b - b / x * x;
    %-----------------%
    
    %-----------------%
    %-detect spindle peak after
    [~, imax(c)] = max(b);
    peak(e, c) = freq.freq(imax(c));
    %-----------------%
    
  end
  %-------------------------------------%
  
end

%-------------------------------------%
%-output and plot
if strcmp(cfg.output, 'average')
  peak = peak(1,:);
end

if cfg.plot

  figure;
  plot(freq.freq, freq.powspctrm)
  legend(data.label{:})
  xlabel('frequency (Hz)')
  ylabel('power spectrum')
  
  hold on
  for c = 1:numel(data.label)
    plot(peak(1,c), freq.powspctrm(c, imax(c)), 'ro')
  end
  
end
%-------------------------------------%
%---------------------------------------------------------%
