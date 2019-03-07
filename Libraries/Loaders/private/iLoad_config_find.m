function loader = iLoad_config_find(loader)
% iLoad_config_find: searches in the config for given names

  if iscellstr(loader)
    loaders = {};
    for index=1:numel(loader)
      this    = iLoad_config_find(loader{index});
      loaders = [ loaders ; this(:) ];
    end
    loader = loaders;
  elseif ischar(loader)
    % test if loader is the user name of a function
    config  = iLoad('','load config');
    formats = config.loaders;
    loaders ={};
    loaders_count=0;
    for index=1:length(formats)
      this_loader = formats{index};
      i1 = strfind(lower(this_loader.name), lower(loader));   if isempty(i1), i1=0; end
      i2 = strfind(lower(this_loader.method), lower(loader)); if isempty(i2), i2=0; end
      i3 = strfind(lower(this_loader.extension), lower(loader)); if iscell(i3), i3 = ~cellfun('isempty', i3); end
      i4 = strfind(lower(this_loader.postprocess), lower(loader)); if iscell(i4), i4 = ~cellfun('isempty', i4); end
      if all(isempty(i3)), i3=0; end
      if all(isempty(i4)), i4=0; end
      if any([ i1 i2 i3 i4 ])
        loaders_count          = loaders_count+1;
        loaders{loaders_count} = this_loader;
      end
    end
    if ~isempty(loaders) loader = loaders; end
  end
  
end % iLoad_config_find

