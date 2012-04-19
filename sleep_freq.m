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