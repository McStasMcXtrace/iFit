function [pars_out,criteria,message,output] = iFunc_feval_list
% get the list of available optimizers and models
  
  output     = {};
  pars_out   = {};
  warn       = warning('off','MATLAB:dispatcher:InexactCaseMatch');
  d = fileparts(which('fminpso'));
  if nargout == 0
    fprintf(1, '\n%s\n', version(iData));
    
    fprintf(1, '      OPTIMIZER DESCRIPTION [%s]\n', [ 'iFit/Optimizers in ' d ]);
    fprintf(1, '-----------------------------------------------------------------\n'); 
  end
  d = dir(d);
  for index=1:length(d)
    this = d(index);
    try
      [dummy, method] = fileparts(this.name);
      options = feval(method,'defaults');
      if isstruct(options) && isfield(options, 'TolFun')
        output{end+1} = options;
        pars_out{end+1}   = method;
        if nargout == 0
          fprintf(1, '%15s %s\n', options.optimizer, options.algorithm);
        end
      end
    end
  end % for
  message  = {}; 
  models   = []; 
  if nargout ~= 1
      % return the list of all available fit functions/models
      d = fileparts(which('gauss'));
      if nargout == 0
        fprintf(1, '\n');
        fprintf(1, '       MODEL DESCRIPTION [%s]\n', [ 'iFit/Models in ' d ]);
        fprintf(1, '-----------------------------------------------------------------\n'); 
      end

      % also search in Specialized and Factory directories
      D = { d, fullfile(d,'Specialized'), fullfile(d,'Factory'), pwd };

      for f_dir = D

        d = dir(f_dir{1});
        for index=1:length(d)
          this = d(index);
          try
            [dummy, method, ext] = fileparts(this.name);
            if strcmp(ext, '.m')
              [mess, options] = evalc([ method '(''identify'')' ]);
            else
              options = [];
            end
            if isa(options, 'iFunc')
              message    = [ message method ];
              models     = [ models options ]; 
              if nargout == 0
                fprintf(1, '%15s %s\n', method, options.Name);
              end
            end
          end
        end % for index(d)
      end % f_dir (model location)
      % sort models with their Name
      [name,index]=sort(get(models,'Name'));
      message  = message(index);
      models   = models(index);
  else
    message  = 'Optimizers and fit functions list'; 
  end
  if nargout == 0 && length(models)
    fprintf(1, '\n');
    % plot all functions
    subplot(models);
  end
  criteria = models;  % return the instantiated models
  warning(warn);


