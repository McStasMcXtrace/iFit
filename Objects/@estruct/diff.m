function structure = diff(A, varargin)
% DIFF find difference between structures
%
% Res=diff(A,B)
%   The result has fields which are either in A and B but not both (difference)
  
  structure = cat(A,varargin{:},'xor');
  
