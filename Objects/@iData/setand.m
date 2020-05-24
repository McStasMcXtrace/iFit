function structure = setand(A, varargin)
% SETAND Merge structures and retain only common fields/properties.
%   Res=SETAND(A,B) The result has fields/properties which are both in
%   A and B (intersection, AND). Object base properties are not
%   affected.
%
% Example: a=iData('a',1:10,'b','blah'); b=iData('a',1,'c',2:50); ...
%   c=setand(a,b); isfield(c, 'a') && ~isfield(c,'b')
% Version: $Date$ $Version$ $Author$
% See also: setcat, setor, setxor

  structure = setor(A,varargin{:},'and');
