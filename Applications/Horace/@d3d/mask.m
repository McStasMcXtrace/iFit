function wout = mask (win, varargin)
% Remove the bins indicated by the mask array
%
% Syntax:
%   >> wout = mask (win, mask_array)
%
% Input:
% ------
%   win                 Input sqw object
%
%   mask_array          Array of 1 or 0 (or true or false) that indicate
%                      which points to retain. 
%                       Numeric or logical array of same number of elements
%                      as the data.
%                       Note: mask will be applied to the stored data array
%                      according as the projection axes, not the display axes.
%                      Thus permuting the display axes does not alter the
%                      effect of masking the data.
%
% Output:
% -------
%   wout                Output dataset.


% Original author: T.G.Perring
%
% $Revision: 301 $ ($Date: 2009-11-03 20:52:59 +0000 (Tue, 03 Nov 2009) $)


% ----- The following shoudld be independent of d0d, d1d,...d4d ------------
% Work via sqw class type

if nargin==1
    wout = win; % trivial case of no sectioning being required
else
    wout = dnd(mask(sqw(win),varargin{:}));
end
