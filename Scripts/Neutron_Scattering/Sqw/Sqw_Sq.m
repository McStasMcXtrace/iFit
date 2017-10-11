function s = Sqw_Sq(s)
% Sqw_Sq: compute the structure factor
%  The structure factor is the integral of the dynamic structure factor along
%  the energy axis. It is representative of the material structure.
%  Its Fourier transform is the pair distribution function g(r).
%
%  The S(q,w) is a dynamic structure factor aka scattering function.
%  This function is basically a call to trapz(s,1)
%
% input:
%   s: Sqw data set e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%
% References: Fischer, Barnes and Salmon, Rev Prog Phys 69 (2006) 233
%
% Example: s = Sqw_Sq(s);
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end
  if ~isa(s, 'iData'), s=iData(s); end

  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      s(index) = feval(mfilename, s(index));
    end
    return
  end
  
  s = Sqw_check(s);
  if isempty(s), return; end

  % compute integral
  s = trapz(s,1); % on q
  
  % reset axes
  title(s,'S(q)');
  if isempty(s.Label), s.Label='S(q)'; end
  

