function [varargout] = verbosity(level, varargin)
% VERBOSITY control the verbosity level of output messages
%
% Set the current verbosity level with e.g.:
%
%  verbosity(level)
%
% Get the current verbosity level:
%
%  level = verbosity; % return the verbosity level (number)
%
% Output a message when it is allowed by the current verbosity level:
%
%  verbosity('normal',  'message')
%  verbosity('verbose', 'message %g %s', arg1, arg2)
%  verbosity(2,         'message %g %s', arg1, arg2) % same as 'info'
%
% Verbosity levels:
%
%  0 error    (no message except errors)
%  1 warning  (important messages not being errors)
%  2 info     (normal messages)
%  3 verbose  (messages for e.g. debugging)
%
% (c) E.Farhi, ILL. License: BSD.

persistent verbosity_level  % stored as a number (faster to process)
persistent verbosity_levels % a list of meaningfull verbosity names

if isempty(verbosity_level) || ~isnumeric(verbosity_level)
  verbosity_level = 2; 
end
if isempty(verbosity_levels)
  verbosity_levels = {'error','warning','info','verbose'}; 
end

if     verbosity_level < 0, verbosity_level=0;
elseif verbosity_level >= numel(verbosity_levels), verbosity_level = numel(verbosity_levels)-1;
end

% verbosity: inquire
if nargin == 0
  % output current verbosity level
  if nargout, varargout(1) = { verbosity_level }; 
  else
    disp([ mfilename ': currently set to ''' verbosity_levels{verbosity_level+1} ''' [' num2str(verbosity_level) ']' ]);
  end
  return
end

% verbosity(level): set level
if nargin >= 1 && ~isnumeric(level)
  switch level
  case {'0','quiet','silent','mute','error','fatal'}
    level = 0;
  case {'1','essential','low','basic','warning','off'}
    level = 1;
  case {'2','info','normal','medium','on'}
    level = 2;
  case {'3','verbose','detail','debug','all','diagnostics','high'}
    level = 3;
  otherwise
    level = 1;
  end
end
if nargin == 1
  verbosity_level = level;
  return
end

% when warning off, supress all messages except 'silent/error' ones. 
q  = warning('query');
vl = verbosity_level;
if numel(q) == 1 && strcmp(q.state, 'off')
  vl = 0;
end

% verbosity(level, message, ...): output
if level <= vl
  v = sprintf(varargin{:});
  if level == 0 % error
    error(v); % so that the error stack is displayed
  end
  disp(v);
end

