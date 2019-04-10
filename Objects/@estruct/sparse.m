function a = sparse(a)
% b = sparse(s) : Convert estruct object storage to sparse
%
%   @estruct/sparse function to use sparse storage, which only stores
%   non zeros elements in Signal, Error and Monitor. This may be usefull
%   for event based storage where most events are zeros. Use estruct/full to
%   revert to full matrix storage.
%
% input:  s: object or array (estruct)
% output: b: object or array (estruct)
% ex:     b=sparse(a);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/full, estruct/pack

a = unary(a, 'sparse');

