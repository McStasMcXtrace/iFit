function a = history(a, meth, varargin)
% HISTORY object command catenation
%   appends the commands 'meth' to the history of 'a' (object or array of)
%   'meth' may be given as a method name plus additional arguments which
%   are transformed into chars to build a command line meth(varargin)

% This file is used by: all methods that record entries in the Command property
%   unary binary axescheck commandhistory copyobj estruct interp setaxis setor subsasgn

% EF 23/09/07 iData implementation

if ~ischar(meth)
  disp([ class(a) '/' mfilename ': command to add in the history should be a char, now ' class(meth) ]);
  return
end

if nargin >= 3 || ~isempty(varargin)
  toadd = ''; tocat = '';
  for i1=1:length(varargin)
    if i1 > 1, c = ','; else c=''; end
    b = varargin{i1};
    if ~isa(b,'estruct') && isstruct(b)
      b = class2str('',b,'eval short');
    end
    if ischar(b)
      if numel(b) > 100, b=[ b(1:20) '...' b((end-20):end) ]; end
      toadd = [ toadd c ' ''' b '''' ];
    elseif isobject(b) && isfield(b, 'Tag')
      toadd = [ toadd c b.Tag ];
      if isempty(tocat), tocat = ' %'; end
      tocat = [ tocat ' <' class(b) ' ' b.Tag ' ' b.Source '> ' ];
    elseif isnumeric(b) || islogical(b)
      if ndims(b) > 2,   b=b(:); end
      if numel(b) > 50, toadd = [ toadd c ' [' sprintf('%g ',double(b(1:20))) '...' sprintf('%g ',double(b((end-20):end))) ']' ];
      elseif isempty(b)
        toadd = [ toadd c ' []' ];
      else
        toadd = [ toadd c ' [' sprintf('%g ',double(full(b))) ']' ];
      end
    else
      toadd = [ toadd c ' <' class(b) '>'  ];
    end
  end
  meth = [ a.Tag '=' meth '(' toadd ');' tocat ];

end

for index=1:numel(a)
  d=a(index);
  if isempty(d.Command),
    d.Command = { meth };
  else
    if ~iscellstr(d.Command), d.Command = { char(d.Command) }; end
    d.Command{end+1} = meth;
  end
  d.Command=d.Command(:);
  % a(index) = d; % not needed as we use handles
end
