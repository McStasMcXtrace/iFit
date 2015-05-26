function s=Sqw_symmetrize(s)
% Sqw_symmetrize(s): extend the S(|q|,w) in both energy sides
%  The resulting S(q,w) is the combination of S(q,w) and S(q,-w), which
%  is thus symmetric in energy:
%     S(q,w) = S(q,-w)
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
%   s:  Sqw data set (classical)
%        e.g. 2D data set with w as 1st axis (rows), q as 2nd axis.
% output:
%   s:  S(|q|,w) symmetrised in energy
%
% Example: Sqw_symmetrize(s, 300)
%
% See also: Sqw_Bosify, deBosify, Sqw_dynamic_range, Sqw_total_xs

  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index)) ];
    end
    s = sqw;
    return
  end

  % check if the data set is Sqw (2D)
  w_present=0;
  q_present=0;
  if isa(s, 'iData') && ndims(s) == 2
    for index=1:2
      lab = lower(label(s,index));
      if any(strfind(lab, 'wavevector')) || any(strfind(lab, 'q')) || any(strfind(lab, 'Angs'))
        q_present=index;
      elseif any(strfind(lab, 'energy')) || any(strfind(lab, 'w')) || any(strfind(lab, 'meV'))
        w_present=index;
      end
    end
  end
  if ~w_present || ~q_present
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be an isotropic S(|q|,w) 2D object. Ignoring.' ]);
    return
  end

  % test if classical
  if ~isempty(findfield(s, 'classical'))
    if s.classical == 0
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be classical. It may already contain the Bose factor in which case the symmetrisation may be wrong.' ]);
    end
  end

  % check if we need to transpose the S(q,w)
  if w_present==2 && q_present==1
    s = transpose(s);
  end

  % test if the data set has single energy side: much faster to symmetrise
  w = s{1}; % should be a row vector
  if isvector(w) && (all(w(:) >= 0) || all(w(:) <= 0))
    signal = get(s, 'Signal');
    signal=[ signal ; signal ];
    [w,index]=unique([ w ; -w ]);
    s{1}=w;
    s = set(s, 'Signal', signal(index,:));
    return
  end
  
  % create a new object with an opposite energy axis
  s_opp = setaxis(s, 1, -s{1});

  % final object (and merge common area)
  s     = combine(s, s_opp);
