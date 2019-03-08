function structure = intersection(A, varargin)
% struct/intersection merge structures and retain only common part
%
% Res=intersection(A,B)
%   The result has fields which are both in A and B (intersection)
  
  structure = cat(A,varargin{:},'and');
  
