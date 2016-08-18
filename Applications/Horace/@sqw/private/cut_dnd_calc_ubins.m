% Create bin boundaries for integration and plot axes from requested limits and step sizes
% Uses knowledge of the range of the data and energy bins of the data to set values for those
% not provided.
%
% Syntax:
%   >> [iiax, ipax, p, noffset, nkeep, mess] = cut_dnd_calc_ubins (pbin, pin, nbin)
%
% Input:
% ------
%   pbin        Cell array of requested limits and bins for plot axes:
%       pbin{1}     Binning along first plot axis
%       pbin{2}     Binning along second plot axis
%                           :
%                   for as many axes as there are plot axes. For each binning entry:
%               - [] or ''          Plot axis: use bin boundaries of input data
%               - [pstep]           Plot axis: sets step size; plot limits taken from extent of the data
%                                   If pstep=0 then use current bin size and synchronise
%                                  the output bin boundaries with the current boundaries. The overall range is
%                                  chosen to ensure that the range of the input data is contained within
%                                  the bin boundaries.
%               - [plo, phi]        Integration axis: range of integration - those bin centres that lie inside this range 
%                                  are included.
%               - [plo, pstep, phi] Plot axis: minimum and maximum bin centres and step size
%                                   If pstep=0 then use current bin size and synchronise
%                                  the output bin boundaries with the current boundaries. The overall range is
%                                  chosen to ensure that the range plo to phi is contained within
%                                  the bin boundaries.
%   pin         Cell array current bin boundaries on plot axes.
%
%   nbin        [2 x ndim] array containing the lower and upper indices of the 
%               extremal elements along each axis
%
% Output:
% -------
%   iax         Index of integration axes into the incoming plot axes  [row vector]
%   iint        Integration range for the new integration axes [(2 x length(iax)) vector]
%   pax         Index of plot axes into the incoming plot axes  [row vector]
%   p           Call array containing bin boundaries along the remaining plot axes [column vectors]
%                   i.e. data.p{1}, data.p{2} ... (for as many plot axes as given by length of data.pax)
%   noffset     Offset along the remaining plot axes of the section from the input signal array (nkeep, below)
%              once the axes to be integrated over have been summed [row vector, length=no. elements of p]
%   nkeep       Section of the input signal array to be retained: [2xndim], ndim=length(pin)
%   mess        Error message; empty if all was OK; non-empty otherwise (in which case all other output are empty)
%               Use this as the sole criterion of succesful operation because empty arrays for other outputs
%               can be valid e.g. iax=[] and iint=[] if all axes are plot axes.
%
% Notes:
%   - if the range of the data is zero along a plot axis, then the two bin boundaries
%    will both have the same value, even if a non-zero bin size was input.
%
%