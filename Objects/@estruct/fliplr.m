function a = fliplr(a)
% b = fliplr(s) : Flip object in left/right direction
%
%   @estruct/fliplr function to flip object in left/right direction
%     With 2D data sets, the X axis (horizontal) is inverted.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=fliplr(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/fliplr, fliplr, estruct/flipud, flipud

a = unary(a, 'fliplr');

