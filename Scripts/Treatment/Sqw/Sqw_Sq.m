function s = Sqw_Sq(s)
% Sqw_Sq: compute the structure factor
%  The structure factor is the integral of the dynamic structure factor along
%  the energy axis. It is representative of the material structure.
%  Its Fourier transform is the pair distribution function g(r).
%
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
% input:
%   s: Sqw data set (non classical, including T Bose factor e.g from experiment)
%        e.g. 2D data set with w as 1st axis (rows), q as 2nd axis.
%
% References: Fischer, Barnes and Salmon, Rev Prog Phys 69 (2006) 233
%
% Example: s = Sqw_Sq(s);


  if nargin == 1, T = []; end

  % handle array of objects
  if numel(s) > 1
    g = [];
    for index=1:numel(s)
      g = [ g feval(mfilename, s(index)) ];
    end
    s = g;
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
    if s.classical == 1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to be classical. The S(q) computation may be wrong. Apply Sqw_Bosify first.' ]);
    end
  end
  
  % check if we need to transpose the S(q,w)
  if w_present==2 && q_present==1
    s = transpose(s);
  end

  % this is the weighting for valid data
  s(s <= 0) = 0; s(~isfinite(s)) = 0;

  % compute integral
  s = trapz(s,1); % on q
  
  % reset axes
  title(s,'S(q)');
  if isempty(s.Label), s.Label='S(q)'; end
  

