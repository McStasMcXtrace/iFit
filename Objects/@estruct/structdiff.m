function structure = structdiff(A, varargin)
% STRUCTDIFF find difference between structures
%
% Res=diff(A,B)
%   The result has fields which are either in A and B but not both (difference, XOR)
%
% Version: $Date$ $Version$ $Author$
% See also: structcat, structintersection, structdiff
  
  structure = structcat(A,varargin{:},'xor');
  
