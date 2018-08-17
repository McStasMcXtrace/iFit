function pars_struct = varargin2struct(pars_names, varargin, flag_lower)
  % varargin2struct: convert a varargin cell into a structure with names parameters
  %
  % pars_struct = varargin2struct(pars_names, varargin, flag_lower)
  % pars_struct = varargin2struct(varargin)
  % pars_struct = varargin2struct(struct)
  %
  % input:
  %   pars_names: expected 'regular' parameters, e.g. {'q', 'T', 'sigma', 'm', 'n', 'dw'}
  %   varargin:   a cell containing input parameters to the function (single cell)
  %   flag_lower: when present and true, all parameters are assumed to be lower case
  % output:
  %   pars_struct: a named structure with parameter values
  %
  % Example:
  %   varargin2struct({'q', 'T'}, {1,2,'dw',5,'q',7})
  %   returns struct(q: 7, T: 2, dw: 5)

  % pars_names = {'q', 'T', 'sigma', 'm', 'n', 'dw'};
  
  pars_struct = struct();
  if nargin == 1
    % [varargin]
    varargin   = pars_names;
    pars_names = {};
    flag_lower = false;
  elseif nargin == 2
    % [names, varargin] or [varargin,flag]
    if iscellstr(pars_names), flag_lower = false;
    else
      flag_lower = varargin;
      varargin   = pars_names;
      pars_names = {};
    end
  end
  
  if flag_lower
    pars_names = lower(pars_names);
  end
  
  % handle case with a struct as input
  if numel(varargin) == 1 && isstruct(varargin{1})
    prop = {};
    reg  = {};
    this = varargin{1};
    for f = fieldnames(this)'
      prop{end+1} = f{1};
      prop{end+1} = this.(f{1});
    end
  else
    [reg, prop] = parseparams(varargin);
  end
  
  % initialize required fields to empty
  for index=1:numel(pars_names)
    pars_struct.(pars_names{index})   = []; % initialize arguments to []
  end
  % transfer 'regular' parameters
  for index=1:numel(reg)
    if numel(pars_names) >= index
      if numel(reg) >= index
        pars_struct.(pars_names{index})   = reg{index};
      end
    end
  end

  % transfer name/value pairs
  for index=1:2:numel(prop)
    if flag_lower
      pars_struct.(lower(prop{index})) = prop{index+1};
    else 
      pars_struct.(prop{index}) = prop{index+1};
    end
  end

