function structure = structintersection(A, varargin)
% structintersection merge structures and retain only common part
%
% Res=structintersection(A,B)
%   The result has fields which are both in A and B (intersection)
  
  structure = structcat(A,varargin{:},'and');
  
