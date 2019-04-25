function structure = structunion(A, varargin)
% STRUCTUNION Concatenate/merge structures fields.
%   STRUCTUNION is equivalent to STRUCTCAT (union, OR).
%
% Version: $Date$ $Version$ $Author$
% See also: structcat, structintersection, structdiff

  structure = structcat(A,varargin{:});