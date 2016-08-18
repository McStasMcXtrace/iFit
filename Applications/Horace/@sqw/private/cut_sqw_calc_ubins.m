% Create bin boundaries for integration and plot axes from requested limits and step sizes
% Uses knowledge of the range of the data and energy bins of the data to set values for those
% not provided.
%
% Syntax:
%   >> [iax, iint, pax, p, urange, mess] = cut_sqw_calc_ubins (urange_in, rot, trans, pbin, pin)
%
% Input:
% ------
%   urange_in   [2x4] array of range of data along the input projection axes (elements must all be finite)
%   rot         Matrix [3x3]     --|  that relate a vector expressed in the
%   trans       Translation [3x1]--|  frame of the bin boundaries to those of urange:
%                                         r'(i) = A(i,j)(r(j) - trans(j))
%   pbin        Cell array of requested limits and bins for integration and plot axes:
%       pbin{1}     Binning along first Q axis
%       pbin{2}     Binning along second Q axis
%       pbin{3}     Binning along third Q axis
%               - [] or ''          Use default bins (bin size and limits)
%               - [pstep]           Plot axis: sets step size; plot limits taken from extent of the data
%               - [plo, phi]        Integration axis: range of integration
%               - [plo, pstep, phi] Plot axis: minimum and maximum bin centres and step size
%
%       pbin{4}     Binning along the energy axis:
%               - omit              Equivalent to [0] and [-Inf,0,Inf]
%               - [] or ''          Use default bins (bin size and limits)
%               - [pstep]           Plot axis: sets step size; plot limits taken from extent of the data
%                                  If pstep=0 then use bin size of energy bins in array en (below) and synchronise
%                                  the output bin boundaries with the reference boundaries. The overall range is
%                                  chosen to ensure that the energy range in urange_in is contained within
%                                  the bin boundaries.
%               - [plo, phi]        Integration axis: range of integration
%           	- [plo, pstep, phi]	Plot axis: minimum and maximum bin centres and step size;
%                                  If pstep=0 then use bin size of energy bins in array en (below) and align
%                                  the output bin boundaries with the reference boundaries. The overall range is
%                                  chosen to ensure that the energy range plo to phi is contained within
%                                  the bin boundaries.
%   pin         Cell array, length 4, of default bin boundaries for each axis. Boundaries assumed equally spaced
%              assumed to be column vectors. 
%               If length(pin{i})==2, will be interpreted as integration range.
%
%   en          Energy bin information used if energy step is zero (see above)
%
%
% Output:
% -------
%   iax         Index of integration axes into the projection axes  [row vector]
%                   e.g. if data is 2D, data.iax=[1,3] means summation has been performed along u1 and u3 axes
%   iint        Integration range along each of the integration axes. [iint(2,length(iax))]
%                   e.g. in 2D case above, is the matrix vector [u1_lo, u3_lo; u1_hi, u3_hi]
%   pax         Index of plot axes into the projection axes  [row vector]
%                   e.g. if data is 3D, data.pax=[1,3,4] means u1, u3, u4 axes are x,y,z in any plotting
%   p           Call array containing bin boundaries along the plot axes [column vectors]
%                   i.e. data.p{1}, data.p{2} ... (for as many plot axes as given by length of data.pax)
%   urange      Array of limits of data that can possibly contribute to the output data structure in the
%               coordinate frame of the output structure [2x4].
%   mess        Error message; empty if all was OK; non-empty otherwise (in which case all other output are empty)
%               Use this as the sole criterion of succesful operation because empty arrays for other outputs
%               can be valid e.g. iax=[] and iint=[] if all axes are plot axes.
%