function a = free(a, varargin)
% b = free(s, parameters, ...) : sets parameter unlock (free) for further fit using the model 's'
%
%   @iFunc/unmlock unlock model parameters during further fit process.
%     to unlock/free a parameter model, use free(model, parameter)
%
%   To unlock/free a set of parameters, you may use a regular expression as:
%     free(model, regexp(model.Parameters, 'token1|token2|...'))
%
%   free(model, {'Parameter1', 'Parameter2', ...})
%     unlock/free parameter for further fits
%   free(model)
%     display freeed parameters
%   free(model,'all')
%     free/unlock all parameters
%   free(model,'none')
%     fix/lock all parameters
%
% input:  s: object or array (iFunc)
%         parameters: names or index of parameters to unlock/free (char or scalar)
% output: b: object or array (iFunc)
% ex:     b=free(a,'Intensity');
%
% Version: $Date$ $Version$ $Author$
% See also iFunc, iFunc/fits, iFunc/mlock, iFunc/xlim, iFunc/fix, iFunc/unmlock

% calls subsasgn with 'free' for each parameter given

a = munlock(a, varargin{:});

if nargout == 0 && ~isempty(inputname(1)) && isa(a,'iFunc')
  assignin('caller',inputname(1),a);
end

