function structure = union(A, varargin)
% struct/union: concatenate/merge structures
%
% Recursively merges fields and subfields of structures A and B to result structure Res
% Simple recursive algorithm merges fields and subfields of two structures
%   Example:
%   A.field1=1;
%   A.field2.subfield1=1;
%   A.field2.subfield2=2;
% 
%   B.field1=1;
%   B.field2.subfield1=10;
%   B.field2.subfield3=30;
%   B.field3.subfield1=1;
% 
%   C=union(A,B);
%
%  by Igor Kaufman, 02 Dec 2011, BSD
% <http://www.mathworks.com/matlabcentral/fileexchange/34054-merge-structures>
  
  structure = cat(A,varargin{:});
  
