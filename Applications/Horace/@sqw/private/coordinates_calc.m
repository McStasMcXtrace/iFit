% Get average values of one or more coordinates for the pixels in each bin of an sqw object
%
% Syntax:
%   >> xvals = coordinates_calc (w, xlist)
%
% Input:
% ------
%   w       sqw object
%   xlist   One or more valid coordinate names (string or cell array of strings)
%           Valid names are:
%               'd1','d2',...       Display axes (for as many dimensions as sqw object has)
%               'h', 'k', 'l'       r.l.u.
%               'E'                 energy transfer
%               'Q'                 |Q|
%
%           E.g. valid xlist include: 'E'  {'d1','d2'}  {'d3','E','Q','k'}, {'h'}
%
% Output:
% -------
%   ok      true if all OK; false if invalid coordinate name given
%
%   mess    empty if all OK, error message otherwise
%
%   xvals   Cell array of xvalues corresponding to the names in xlist, with size
%          of each entry being that of the signal array.
%           If input was a single name as a character string, then xvals is a numeric array; if
%          a single name in a cell array, then xvals is a cell array with one element.
%
%   xpix    Cell array of corresponding values for each pixel (column vectors)
%
%   xvar    Variance of xvalues, size equal to that of the signal array.
%
%   xdevsqr Cell array of corresponding squared-deviation for each pixel (column vectors)
%