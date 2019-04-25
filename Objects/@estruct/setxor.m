function structure = setxor(A, varargin)
% SETXOR Find difference between structures fields.
%   Res=SETXOR(A,B) The result has fields which are either in A and
%   B but not both (difference, XOR).
%
% Version: $Date$ $Version$ $Author$
% See also: setcat, structintersection, setor, setand

  structure = setor(A,varargin{:},'xor');
