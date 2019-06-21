function t = title(a, varargin)
% TITLE Signal title.
%   TITLE(s) returns the current Signal title. This is equivalent to
%   LABEL(s, 0).
%
%   TITLE(s, 'text') sets  the Signal title. This is equivalent to
%   LABEL(s, 0, 'text').
%
%   The object name is obtained with s.Name
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot

t = label(a, 0, varargin{:});
if isempty(t), t = getalias(a,'Signal'); end
if isempty(t) || ~ischar(t), t = a.Name; end
