function structure = setand(A, varargin)
% SETAND Merge structures and retain only common fields.
%   Res=SETAND(A,B) The result has fields which are both in
%   A and B (intersection, AND).
%
% Version: $Date$ $Version$ $Author$
% See also: setcat, setor, setxor

  structure = setor(A,varargin{:},'and');
