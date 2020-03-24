function loaders = iLoad_config_find(loader, config)
% iLoad_config_find: searches in the config for given names

  if nargin < 2
    config = iLoad('','load config');
  end
  if iscellstr(loader)
    loaders = {};
    for index=1:numel(loader)
      this    = iLoad_config_find(loader{index}, config);
      loaders = [ loaders ; this(:) ];
    end
  elseif ischar(loader)
    % test if loader is the user name of a function

    loaders = {};
    for index=1:length(config.loaders)
      this = config.loaders{index};
      keepme=0;
      if ~isempty(strfind(lower(this.name), lower(loader))) | ~isempty(strfind(lower(this.method), lower(loader)))
        keepme = 1;
      elseif isfield(this,'extension') && any(strncmpi(loader, this.extension, length(loader)))
        keepme = 1;
      end
      if keepme, loaders = { loaders{:} this }; end
    end
    if length(loaders) == 1
      loaders = loaders{1}; 
    end
  end
  
end % iLoad_config_find

