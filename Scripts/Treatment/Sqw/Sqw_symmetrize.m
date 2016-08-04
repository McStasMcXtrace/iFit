function s=Sqw_symmetrize(s)
% Sqw_symmetrize(s): extend the S(|q|,w) in both energy sides
%  The resulting S(q,w) is the combination of S(q,w) and S(q,-w), which
%  is thus symmetric in energy:
%     S(q,w) = S(q,-w)
%
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
%  The incoming data set should NOT contain the Bose factor, that is it
%    should be 'classical'.
%  To obtain a 'classical' S(q,w) from an experiment, use first:
%    Sqw_deBosify(s, T)
%
% The positive energy values in the S(q,w) map correspond to Stokes processes, 
% i.e. material gains energy, and neutrons loose energy when scattered.
%
% input:
%   s:  Sqw data set (classical, often labelled as S*)
%        2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
% output:
%   s:  S(|q|,w) symmetrised in energy
%
% Example: Sqw_symmetrize(s, 300)
%
% See also: Sqw_Bosify, deBosify, Sqw_dynamic_range, Sqw_scatt_xs
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end
  if ~isa(s, 'iData'), s=iData(s); end

  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index)) ];
      s(index) = iData;
    end
    s = sqw;
    return
  end

  s = Sqw_check(s); % in private
  if isempty(s), return; end

  % test if classical
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if s.classical == 0
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be classical. It may already contain the Bose factor in which case the symmetrisation may be wrong.' ]);
    end
  else
    s     = setalias(s,'classical', 1);
  end

  % test if the data set has single energy side: much faster to symmetrise
  w = s{1}; % should be a row vector
  w0= w;
  if isvector(w) && (all(w(:) >= 0) || all(w(:) <= 0))
    signal = get(s, 'Signal');
    signal=[ signal ; signal ];
    [w,index]=unique([ w ; -w ]);
    s{1}=w;
    s = set(s, 'Signal', signal(index,:));
    clear signal
    
    if ~isempty(getalias(s,'Error')) && ~strcmp(getalias(s,'Error'),'sqrt(this.Signal)')
       err = get(s, 'Error');
       if all(err(:) > 0 )
         err=[ err ; err ];
         [w,index]=unique([ w0 ; -w0 ]);
         s = set(s, 'Error', err(index,:));
       end
    end
    clear err

    if ~isempty(getalias(s,'Monitor')) && ~isscalar(get(s, 'Monitor'))
      m = get(s, 'Monitor');
      if all(m(:) > 0 )
        m=[ m ; m ];
        [w,index]=unique([ w0 ; -w0 ]);
        s = set(s, 'Monitor', m(index,:));
      end
    end
    clear m
    
    return
  end
  
  % create a new object with an opposite energy axis

  % final object (and merge common area)
  s     = combine(s, setaxis(s, 1, -s{1}));
  
