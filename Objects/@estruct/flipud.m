function a = flipud(a)
% b = flipud(s) : Flip object in up/down direction
%
%   @estruct/flipud function to flip object in up/down direction
%     With 2D data sets, the Y axis (vertical) is inverted.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=flipud(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/fliplr, fliplr, estruct/flipud, flipud

a = unary(a, 'flipud');

