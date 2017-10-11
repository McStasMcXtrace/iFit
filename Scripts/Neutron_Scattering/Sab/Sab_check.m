function s = Sab_check(s)
% Sab_check: check if a 2D iData is a S(alpha,beta).
%
% conventions:
% w = omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% input:
%   s:  Sqw data set e.g. 2D data set with beta as 1st axis (rows), alpha as 2nd axis (columns).
%
% Example: sqw=Sqw_check('SQW_coh_lGe.nc');
%
% See also: Sqw_Sab, Sqw_check
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end
  
  % handle array of objects
  if numel(s) > 1
    sab = [];
    for index=1:numel(s)
      sab = [ sab feval(mfilename, s(index)) ];
    end
    s(index)=iData; % free memory
    s = sab;
    return
  end
  
  if isempty(s), return; end
  
  % check if the data set is Sqw (2D)
  alpha_present=0;
  beta_present=0;
  q_present=0;
  w_present=0;
  if isa(s, 'iData') && ndims(s) == 2
    for index=1:2
      lab = lower(label(s,index));
      def = getaxis(s, num2str(index));
      if ischar(def), lab = [ def ' ' lab ]; end
      if isempty(lab), lab=lower(getaxis(s, num2str(index))); end
      lab = strread(lab, '%s'); % split string into cell
      if strcmpm(lab, {'alpha','a'}) % strcmpm = multiple strcmpi is private below
        alpha_present=index;
      elseif strcmpm(lab, {'beta','b'})
        beta_present=index;
      elseif strcmpm(lab, {'wavevector','momentum','q','k','angs'})
        q_present=index;
      elseif strcmpm(lab, {'energy','frequency','w','e','mev'})
        w_present=index;
      end
    end
  end
  
  % conversions
  if q_present && w_present && (~alpha_present || ~beta_present)
    s = Sqw_Sab(s);
    s = Sab_check(s);
    return
  end
  
  % search for Sab parameters
  if ~isfield(s, 'parameters')
    s = Sqw_parameters(s, 'Sab');
  end
  if ~alpha_present || ~beta_present
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
    disp('    does not seem to be an isotropic S(alpha,beta) 2D object. Ignoring.');
    disp(' ndims alpha  beta     q     w    ')
    disp([ ndims(s) alpha_present beta_present q_present w_present ]);
    s = [];
    return
  end

  % check if we need to transpose the S(q,w)
  if alpha_present==1 && beta_present==2
    s = transpose(s);
  end
  
  % this is the weighting for valid data
  s(~isfinite(s)) = 0;
  
  % check 'classical' and 'symmetric'
  if isfield(s,'classical') 
    classical0 = s.classical;
  else
    classical0 = [];
  end
  
  % handle ENDF keywords
  if isfield(s,'LASYM')
    s.classical = ~s.LASYM;
    classical0  = s.classical;
  end

% ------------------------------------------------------------------------------
function flag=strcmpm(str, words)
% multiple strcmp
%
% input:
%   str:   string or cellstr of tokens
%   words: string or cellstr of words to search

  flag = false;
  if ischar(str),   str=strread(str, '%s','delimiter',' ,; $()[]{}=|<>&"/\:"'''); end
  if ischar(words), words = strread(words, '%s'); end;
  for index=1:numel(words)
    if any(strcmpi(str, words{index}))
      flag = true; return;
    end
  end
