function structure = setxor(A, varargin)
% SETXOR Find difference between structures fields/properties.
%   Res=SETXOR(A,B) The result has fields/properties which are either in A and
%   B but not both (difference, XOR). Object base properties are not
%   affected.
%
% Example: a=estruct('a',1:10,'b','blah'); b=estruct('a',1,'c',2:50); ...
%   c=setxor(a,b); ~isfield(c, 'a') && isfield(c,'b')
% Version: $Date$ $Version$ $Author$
% See also: setcat, structintersection, setor, setand

  structure = setor(A,varargin{:},'xor');
