function [pks,locs,extra] = findpeaks(a, dim, varargin)
% FINDPEAKS Find local maxima.
%   PKS = FINDPEAKS(A) returns a vector with the local maxima (peaks) of the  
%   input object.
%
%   FINDPEAKS(A, DIM) does the same along dimension DIM with the object 
%   projection. Specifying DIM=0 works on all dimensions (default).
%
%   FINDPEAKS(A,DIM,'PARAM1',VALUE1,...) specifies additional parameters for the 
%   search:
%     MinPeakHeight:   Minimum peak height (positive scalar, 2*std(A))
%     MinPeakDistance: Minimum separation between (positive integer, 4 pts)
%     MinPeakWidth:    Minimum width of peaks (positive integer, 2 pts)
%     DoubleSided:     True if data takes positive and negative values
%                      (boolean, true).
%
%   [PKS,LOCS,EXTRA] = FINDPEAKS(A,...) returns the peak heights PKS, locations 
%   LOCS, and additional information in structure EXTRA.
%     EXTRA.sign:     1 for maximum, -1 for minimum.
%     EXTRA.height:   peak height (signal value at peak), same as PKS.
%     EXTRA.indices:  indices of peak locations in the signal.
%     EXTRA.baseline: baseline below the signal (when noisy, sharp peaks).
%   LOCS is a single vector for 1D objects or when DIM is specified, and a cell 
%   array of axis values for nD objects.
%
% References: Slavic, NIM 112 (1973) 253 ; M. Morhac, NIM A 600 (2009) 478 
%
% Example: t = 2*pi*linspace(0,1,1024)'; ...
%   y = sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + ...
%   0.1*sin(15.3*t+1/3); a=iData(y); a{1}=t; ...
%   [pks x extra] = findpeaks(a); ...
%   plot(a,'hide-errorbars'); hold on; plot(x, pks, 'ro'); ...
%   delete(gcf); numel(pks)==6
% Example: a=iData(peaks); ...
%   [pks x extra] = findpeaks(a); ...
%   plot(a); hold on; plot3(x{2},x{1},ones(size(x{1})),'ro'); ...
%   delete(gcf); numel(pks)>30
% Version: $Date$ $Version$ $Author$
% See also iData, iData/median, iData/mean, iData/std

% inline functions: BaseLine, PeakWidth
% private: findpeaks_octave, from Octave.

  if nargin < 2,   dim=0; end
  if isempty(dim), dim=0; end
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

  if abs(dim) > prod(ndims(a))
    dim = 1;
  end

  % we first compute projection of iData on the selected dimension
  b = copyobj(a); set(b,'Error',0);
  if dim > 0 && ndims(a) > 1 && ~isvector(a)
    b = camproj(b, dim);
  else
    if dim <= 0
      b = subsref(b, substruct('()',':'));
    end
  end
  
  % first we make sure the axis is regularly sampled
  b = unique(b); % unique and sorted axes
  signal    = private_cleannaninf(get(b,'Signal')); 
  % if dim <= 0, signal = signal(:); end
  
  warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
  [pks,locs,extra] = findpeaks_octave(signal, varargin{:});
  
  extra.indices   = locs;
  extra.baseline  = BaseLine(signal, 5);
  extra.width     = PeakWidth(signal, 5); 
  extra.sign      = findpeaks_slavic(signal);
  extra.sign      = extra.sign(locs);

  % compute axes values 
  if dim, c=b; else c=a; end
  x=cell(1,ndims(c));
  [x{:}] = ind2sub(size(c),locs); % [ I, J, K, ..] indices in A signal
  for index=1:ndims(c)
    ax = getaxis(c,index);
    if numel(ax) == size(c,index)
      x{index} = ax(x{index});
    elseif numel(ax) == prod(size(a))
      x{index} = ax(locs);
    end
  end
  locs = x;
  if numel(locs) == 1, locs = locs{1}; end

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
