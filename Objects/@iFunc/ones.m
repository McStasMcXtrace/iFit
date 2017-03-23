function s = ones(iFunc_in,varargin)
% s = ones(s,N,M,P,...) : initialize an iFunc array
%
%   @iFunc/ones function to create an array of 's' iFunc objects
%   The object 's' is duplicated into an array. Use s=iFunc to get an empty array.
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex: ones(iFunc,5,5) will create a 5-by-5 empty iFunc array
%     ones(s,5,5) will return a 5-by-5 array filled with 's'
%
% Version: $Date$
% See also iFunc

% EF 27/07/00 creation
% EF 23/09/07 iFunc impementation

if nargin == 1
    s = zeros(iFunc_in, size(iFunc_in));
else
    s = zeros(iFunc_in, varargin{:});
end

