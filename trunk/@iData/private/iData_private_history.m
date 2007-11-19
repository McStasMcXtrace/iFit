function a = iData_private_history(a, meth, obj, varargin)
% iData_private_history: iData command catenation
%   appends the commands 'meth' to the history of 'a' (iData or cell/array of iData)
%   meth may be given as a method name plus additional arguments which
%   are transformed into chars to build a command line meth(varargin)

% EF 23/09/07 iData impementation

if ~ischar(meth)
  disp('iData/private/iData_private_history: command to add in the history should be a char/cellstr');
  return
end

if nargin >= 3 | length(varargin)
  toadd = '';
  for i1=1:length(varargin)
    b = varargin{i1};
    if ischar(b)
      toadd = [ toadd ', ''' b '''' ];
    elseif isnumeric(b) | islogical(b) 
      if ndims(b) > 2, b=b(:); end
      if numel(b) > 100, b=b(1:100); end 
      toadd = [ toadd ', ' mat2str(double(b)) ];
    elseif isa(b, 'iData')
      toadd = [ toadd ', <' class(b) ' ' b.Tag ' ' b.Source '>'  ];
    else
      toadd = [ toadd ', <' class(b) '>'  ];
    end
  end
  meth = [ a.Tag '=' meth '(' obj.Tag toadd ');' ];
end

meth = cellstr(meth);
for index=1:length(a)
  d=a(index);
  if isempty(d.Command), d.Command = ''; end
  if ~iscellstr(d.Command), d.Command = cellstr(d.Command); end
  d.Command = vertcat(d.Command(:), meth(:));
  a(index) = d;
end


