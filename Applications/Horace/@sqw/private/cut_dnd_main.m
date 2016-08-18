% Take a cut from an dnd-type sqw object by integrating over one or more of the plot axes.
% 
%   >> w = cut_dnd_main (data_source, ndims, p1_bin, p2_bin...)     
%                                           % cut plot axes, keeping existing integration ranges
%                                           % (as many binning arguments as there are plot axes)
%   >> w = cut_dnd_main (..., '-save')      % Save cut to file (prompt for output file)
%   >> w = cut_dnd_main (...,  filename)    % save cut to named file
%
%   >> cut(...)                     % save cut to file without making output to workspace 
% 
% Input:
% ------
%   data_source_in  Data source: dnd object or filename of a file with sqw-type or dnd-type data
%                  (character string or cellarray with one character string)
%
%   ndims           Number of dimensions of the sqw data
%
%   p1_bin          Binning along first plot axis
%   p2_bin          Binning along second plot axis
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
%
% Output:
% -------
%   w              Output data object
%