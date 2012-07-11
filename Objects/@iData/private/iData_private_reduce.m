function a=iData_private_reduce(a, max_length)
% rebinning object so that its number of elements is smaller than 1e6

if nargin == 1
  % how many axes are 'big' ? if more than one, set max_length to 100
  big_axes   = 0;
  max_length = 1000;
  for index=1:ndims(a)
    if length(getaxis(a, index)) > max_length
      big_axes = big_axes +1;
    end
  end
  if big_axes > 1, max_length = 100;
  else             max_length = 1000;
  end
end

S.type='()';
S.subs=cell(1,ndims(a));

% scan dimensions and reduce them to max_length elements per dimension
for index=1:ndims(a)
  lx = length(getaxis(a, index));
  if lx > max_length
    S.subs{index}=ceil(linspace(1,lx,max_length));
  else
    S.subs{index} = ':';
  end
end

% perform rebinning
a = subsref(a, S);

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),a);
end

