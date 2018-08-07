function  T = Sqw_getT(s, prop, raw)
% Sqw_getT: search for a property value in a data set
%
%   T = Sqw_getT(s, prop, raw)
%
%   T = Sqw_getT(s, prop, 'Bose')
%     compute the temperature from the Bose factor (detailed balance).
%
% input:
%   s: any iData object, including S(q,w) and DOS ones.
%   prop: a list of equivalent properties (numeric) to search for.
%   raw: optional, when specified, use raw output, else compute the mean value.
%
% output:
%   T: the property value (mean value when 'raw' not specified)

  if nargin < 2,    prop = []; end
  if nargin ==3 && strcmpi(raw, 'Bose')
    raw = 'Bose';
  elseif nargin < 3
    raw = false; 
  else raw = true; end
  
  T = [];
  if strcmp(raw, 'Bose')
    prop= {'Temperature','T'};
    T   = Sqw_getT_Bose(s);
    if ~isempty(T), return; end
  end
  
  if isempty(prop), prop={'Temperature','T'}; end
  if ~iscellstr(prop), prop = { prop }; end
  
  % handle arrays
  if numel(s) > 1
    for index=1:numel(s)
      t = Sqw_getT(s(index), prop);
      if isempty(t), t=nan; end
      T = [ T t ];
    end
    return
  end
  
  for field=prop
    % does it already exist as a property ?
    if isfield(s,field{1}), T=get(s,field{1}); end
    
    % search case sensitive, the first match only
    if isempty(T) || all(T(:)<0) || ~isnumeric(T) || ~isvector(T)
      f = findfield(s,field{1},'exact numeric first');
      if ~isempty(f), T = get(s,f); end
    end
    
    % not found ? search non case sensitive, all matches
    if isempty(T) || all(T(:)<0) || ~isnumeric(T) || ~isvector(T)
      f = findfield(s,field{1},'numeric cache');
      for index=1:numel(f)
        try
          T = get(s,f{1});
          if ~isempty(T) && isnumeric(T) && isvector(T) && all(T(:)>0), break; end
        end
      end
    end

    if ~isempty(T) && isnumeric(T) && isvector(T) && all(T(:)>0), break; end
  end
  
  if ~raw
    if isvector(T), T = mean(T(:)); else T=[]; end % also for scalars
  end
  
% ------------------------------------------------------------------------------
  
function T = Sqw_getT_Bose(s)
  % Compute the temperature from the Bose factor
  
  w = getaxis(s,1);
  T = [];
  
  if any(w(:) < 0) && any(w(:) > 0)
    % restrict the energy axis to the common +/- range
    w1 = max(w(:)); w2 = max(-w(:)); w_max = min([w1 w2]);

    if w1 ~= w_max || w2 ~= w_max
      s_res  = ylim(s, [-w_max w_max]); % restricted to [-w:w] range
    else
      s_res = s;
    end
    % get axes
    w = getaxis(s_res,1);
    
    % we compare the s(q,w) and s(q,-w)
    s_opp = setaxis(s_res, 1, -w);
    s_opp = sum(s_opp,2); s_opp = sort(s_opp, 1);

    s_res = sum(s_res,2); s_res = sort(s_res, 1);

    % the ratio should be S(q,w)/S(q,-w) = exp(hw/kT)
    % so log(S(q,w)) - log(S(q,-w)) = hw/kT
    wrn = warning;
    warning off;
    log_s_ratio = log(s_res) - log(s_opp);
    warning(wrn);
    w = getaxis(log_s_ratio,1);
    clear s_res s_opp
    
    % mean_log_ratio = mean(log_s_ratio,0);
    % std_log_ratio  = std(log_s_ratio,0);
    
    % compute the temperature from the Data
    % log_s_ratio should be a constant if S(q,w) contains Bose
    % then kT = w./log_s_ratio
    T         = w./log_s_ratio*11.6045; % 1 meV = 11.6045 K
    T         = getaxis(T,0);
    T         = T(isfinite(T));
    index = find(T ~= 0);
    if ~isempty(index)
      T = T(index);
    end
    T         = mean(abs(real(T)));
  end
 
