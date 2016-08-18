% Make constant interval bin boundaries
%
% Required input:
% ---------------
%   pbin    Binning argument: [plo, pstep, phi].
%
%           If pstep>0: defines bin centres and bin size
%               If plo and phi both finite, then bin centre alignment is to plo, and phi is 
%                   interpreted as the data limit to be included in the bins
%               If plo=-Inf, then bin centre alignment is to phi, and lowest bin boundary
%                   set to include lowest data point given by min(range) (below)
%               If phi=+Inf, then bin centre alignment is to plo, and highest bin boundary
%                   set to include highest data point given by max(range) (below)
%               If plo=-Inf and phi=+Inf, then bin centre alignment is to zero, and
%                   lower and upper bin boundaries set to include range of data set by range (below).
%
%           If pstep=0: bin size and alignment will be determined from optional input argument
%               pref (below), and plo, phi will be interpreted as ranges of data to be covered by bins.
%               If either plo=-Inf or phi=+Inf, then corresponding range will be taken from range (below).
%           
% Optional input:
% ---------------
%   range   Range of data to be covered by bins.
%            - can be simply the lower and upper data ranges: [range_lo, range,hi]
%            - more generally, array of data points, assumed to be in increasing order
%           Used to set the overall extent of the bin boundaries when the corresponding lower and/or upper
%          limit(s) in pbin are -Inf and/or Inf. (Argument is ignored if lower and upper limits in pbin are
%          finite, so can set to anything in this case)
%
%   pref    Reference bin size and alignment information used when pstep in pbin is zero.
%           If pref scalar: pstep set to pref, bin centres aligned to zero
%                (i.e. pref is bin size)
%           If pref array:  pstep=pref(2)-pref(1), bin boundaries aligned to pref(1).
%                (i.e. pref can be a set of bin boundaries, and the output will be aligned to pref).
%
%   inside  If pstep=0, so that pref was used to determine the binning size and alignment, then
%          if inside=true, only bins whose centres in the range determined by pbin and range will be retained
%          This inverts the interpretation of the range - useful for energy binning arguments.
%
% Ouput:
% ------
%   p       Column vector of bin boundaries
%