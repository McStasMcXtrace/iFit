function a = uplus(a)
%  +  Unary plus (copyobj).
%    +A for objects is A. A new object is created, with same content,
%    but different Tag/ID and Date (same as COPYOBJ).
%
%    B = UPLUS(A) is called for the syntax '+A'.
%
% Example: a=iData(peaks); ~isempty(+a)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

a = unary(a, 'uplus');
