function data = iLoad_display_formats(data, config)

    
    fprintf(1, ' EXT                    READER  DESCRIPTION [%s]\n', 'iLoad');
    fprintf(1, '-----------------------------------------------------------------\n');
    % scan loaders
    if isstruct(data) && isfield(data,'loaders'), data=data.loaders; end
    for index=1:length(data)
      if iscell(data) this=data{index}; 
      elseif isstruct(data) this=data(index); 
      else continue; end
      % build output strings
      if isfield(this,'postprocess'), 
        if ~isempty(this.postprocess)
          if iscell(this.postprocess)
            for i=1:length(this.postprocess)
              this.method = [ this.method '/' this.postprocess{i} ]; 
            end
          else
            this.method = [ this.method '/' this.postprocess ]; 
          end
        end
      end
      if length(this.method)>25, this.method = [ this.method(1:22) '...' ]; end
      if ~isfield(this,'extension'), this.extension = '*';
      elseif isempty(this.extension), this.extension='*'; end
      % display
      pat = this.patterns;
      if iscellstr(pat) && numel(pat) > 0, pat = sprintf('"%s" [%i]', pat{1},numel(pat));
      elseif ~ischar(pat), pat='';
      elseif isempty(pat), pat='[text]';
      else pat = sprintf('"%s"', pat); end
      if iscellstr(this.extension)
        fprintf(1,'%4s %25s  %s %s\n', upper(this.extension{1}), this.method,this.name, pat);
        for j=2:length(this.extension),fprintf(1,'  |.%s\n', upper(this.extension{j})); end
      else
        fprintf(1,'%4s %25s  %s %s\n', upper(this.extension), this.method,this.name,pat);
      end

    end
    disp([ '% iLoad configuration file: ' config.FileName ]);
    data = data(:)';
