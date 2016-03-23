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

  if nargin == 0, return; end
  if ~isa(s, 'iData')
    disp([ mfilename ': ERROR: The data set should be an iData object, and not a ' class(s) ]);
    return; 
  end
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
  
  s = Sqw_check(s);
  if isempty(s), return; end

  % test if classical
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if s.classical == 1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to be classical. The S(q) computation may be wrong. Apply Sqw_Bosify first.' ]);
    end
  end

  % compute integral
  s = trapz(s,1); % on q
  
  % reset axes
  title(s,'S(q)');
  if isempty(s.Label), s.Label='S(q)'; end
  

