function a = uplus(a)
%  +  Unary plus.
%    +A for objects is A. A new object is created, with same content,
%    but different Tag/ID and Date (same as COPYOBJ).
%
%    B = UPLUS(A) is called for the syntax '+A'.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/uminus, estruct/abs, estruct/real, estruct/imag, estruct/uplus

a = unary(a, 'uplus');