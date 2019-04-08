function structure = structintersection(A, varargin)
% STRUCTINTERSECTION merge structures and retain only common part
%
% Res=structintersection(A,B)
%   The result has fields which are both in A and B (intersection, AND)
%
% Version: $Date$ $Version$ $Author$
% See also: structcat, structintersection, structdiff
  
  structure = structcat(A,varargin{:},'and');
  
