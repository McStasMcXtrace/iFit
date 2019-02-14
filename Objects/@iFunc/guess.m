function [p, model, ax, name] = guess(s, varargin)
% b = guess(s) : guess starting parameters for model
%
%   @iFunc/guess makes a guess of the starting parameters for model
%     using default or optional  axes and signal
%
%  [p, model, ax, name] = guess(s)
%     returns as well the updated model and guesses axes.
%
%  [p, model, ax, name] = guess(s, x,y,..., signal)
%     specifies optional input axes and signal.
%
% syntax:
%   p=guess(s)
%   p=guess(s, x,y,...)
%   p=guess(s, x,y,..., signal)
%
% input:  s: object or array (iFunc)
% output: p: parameter values
% ex:     p=guess(gauss);
%
% Version: $Date$
% See also iFunc, iFunc/feval

p = feval(s, 'guess_only', varargin{:});

if ~isempty(inputname(1))
  assignin('caller',inputname(1),s);
end
