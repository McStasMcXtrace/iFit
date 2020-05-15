function a = uplus(a)
%  +  Unary plus (copyobj).
%    +A for objects is A. A new object is created, with same content,
%    but different Tag/ID and Date (same as COPYOBJ).
%
%    B = UPLUS(A) is called for the syntax '+A'.
%
% Example: a=estruct(peaks); ~isempty(+a)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

a = unary(a, 'uplus');
