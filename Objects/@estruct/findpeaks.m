function [pks,locs,w,extra] = findpeaks(a, dim, varargin)
% FINDPEAKS Find local maxima.
%   PKS = FINDPEAKS(A) returns a vector with the local maxima (peaks) of the  
%   input object along first dimension.
%
%   FINDPEAKS(A, DIM) does the same along dimension DIM.
%
%   FINDPEAKS(A, DIM, 'PARAM1',VALUE1,...) specifies additional parameters for the 
%   search:
%     MinPeakHeight:   Minimum peak height (positive scalar, 2*std(A))
%     MinPeakDistance: Minimum separation between (positive integer, 4 pts)
%     MinPeakWidth:    Minimum width of peaks (positive integer, 2 pts)
%     DoubleSided:     True if data takes positive and negative values
%                      (boolean, true).
%
%   [PKS,LOCS,W,EXTRA] = FINDPEAKS(A,...) returns the peak heights PKS, locations 
%   LOCS, peak widths W and additional information in structure EXTRA.
%     EXTRA.sign:     1 for maximum, -1 for minimum.
%     EXTRA.height:   peak height (signal value at peak), same as PKS.
%     EXTRA.width:    peak width, same as W.
%     EXTRA.indices:  indices of peak locations in the signal.
%     EXTRA.baseline: baseline below the signal (when noisy, sharp peaks).
%
% References: Slavic, NIM 112 (1973) 253 ; M. Morhac, NIM A 600 (2009) 478 
%
% Example: t = 2*pi*linspace(0,1,1024)'; ...
%   y = sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + ...
%   0.1*sin(15.3*t+1/3); a=estruct(t,y); ...
%   [pks x w extra] = findpeaks(a); ...
%   plot(a); hold on; plot(x, pks, 'ro'); delete(gcf); numel(pks)==6
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/median, estruct/mean, estruct/std

% inline functions: BaseLine, PeakWidth
% private: findpeaks_octave, from Octave.

  if nargin < 2, dim=1; end
  if numel(a) > 1
    pks = cell(size(a)); locs = pks; w = pks; extra = pks;
    for index=1:numel(a)
      [p,l,v,e] = peaks(a(index), dim, varargin{:});
      pks{index}   = p;
      locs{index}  = l;
      w{index}     = v;
      extra{index} = e;
    end
    return
  end

  if dim == 0 || abs(dim) > prod(ndims(a))
    dim = 1;
  end

  % we first compute projection of estruct on the selected dimension
  if ndims(a) > 1 && ~isvector(a)
    b = camproj(a, abs(dim));
  else
    b = copyobj(a);
  end
  
  % first we make sure the axis is regularly sampled
  b = unique(b); % unique and sorted axes
  signal    = private_cleannaninf(get(b,'Signal')); 
  x         = getaxis(b,1);
  
  warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
  [pks,locs,extra] = findpeaks_octave(signal, varargin{:});
  
  extra.indices   = locs;
  extra.baseline  = BaseLine(signal, 3);
  extra.width     = PeakWidth(signal-extra.baseline, 3)*min(diff(x));
  extra.sign      = findpeaks_slavic(signal);
  extra.sign      = extra.sign(locs);
  extra.width     = extra.width(locs);
  w               = extra.width;
  locs            = x(locs); % now as X values, not indices.
end

% ==============================================================================
% inline functions: BaseLine, PeakWidth
% ==============================================================================
function index=findpeaks_slavic(signal)
% Slavic, NIM 112 (1973) 253 ; M. Morhac, NIM A 600 (2009) 478 
  Gmm = circshift(signal, -1); Gmm((end-1):end) = 0;
  Gpm = circshift(signal,  1); Gpm(1) = 0; 
  % a Max is such that Gmm < signal & Gpm < signal
  index_max = find(Gmm < signal & Gpm < signal);
  % a Min is such that Gmm > signal & Gpm > signal
  index_min = find(Gmm > signal & Gpm > signal);
  index = zeros(size(signal));
  index(index_max)= 1;
  index(index_min)=-1;
end

function baseline = BaseLine(y, m)
% BaseLine: compute signal baseline from: M. Morhac, NIM A 600 (2009) 478
% computes baseline estimate along signal 'y' of given length 'n'
% with contrast parameter 'm' (see paper above).

  baseline = [];
  if nargin == 1, m=0; end
  if (m<=0) m=ceil(max(5, length(y)/50)); end % automatic largest width estimate
  if (length(y)<=m) return; end
  
  % improved from 'Algorithm C' from M. Morhac, NIM A 600 (2009) 478
  baseline=y;
  shifts = zeros(m+1, length(y)-m);
  new_length=length(y)-m;
  for i=0:m;
    shifts(i+1, :) = y((i+1):(i+new_length));
  end
  shifts = min(shifts);
  i = ceil(m/2);
  baseline((i+1):(i+new_length)) = shifts;
end

function sigma = PeakWidth(signal, m) 
% PeakWidth: estimate peak width from: M. Morhac, NIM A 600 (2009) 478
% computes peak width estimate sigma along signal of given length
% with contrast parameter 'm' (see paper above).
% the peak width is given in bins.

  sigma = [];
  sz = size(signal);
  signal = signal(:);
  if nargin == 1, m=0; end
  
  if (m<1) m=ceil(max(5, length(signal)/50)); end
  if length(signal) < 4*m, m=3; end
  if (length(signal)<=m) return; end
  
  % Gaussian product function as of Slavic, NIM 112 (1973) 253 
  Gmm = circshift(signal, -m); Gmm((end-m):end) = 0;
  Gpm = circshift(signal,  m); Gpm(1:m) = 0; 
  
  sigma = ones(sz)*max(signal);
  index = find(Gmm < signal & Gpm < signal & Gmm~=0 & Gpm~=0);
  sigma(index) = signal(index).*signal(index)./Gmm(index)./Gpm(index);
  index= find(sigma>1);
  sigma(index) = m./sqrt(log(sigma(index)));
  sigma = reshape(sigma, sz);
end
