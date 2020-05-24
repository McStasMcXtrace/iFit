function vout = arrayfun(fun, s)
% ARRAYFUN Apply a function to each element of iData array.
%    A = ARRAYFUN(FUN, B) applies the function specified by FUN to each
%    element in iData array B. To use a class method explicitly, specify e.g.
%      FUN = @B.method or @B.method(arg1, ...)
%
% Example: s = iData(1:10); numel(arrayfun(@(s)size(s,2),[s s])) == 2
% Version: $Date$ $Version$ $Author$
% See also cellfun, structfun, function_handle, cell2mat

% get the number of output arguments from 'fun'
  vout        = cell(1,numel(s));
  nout        = nargout(fun);

% call FUN for each array element
  for index=1:numel(s)
    if nout>0
      this_out=cell(1,nout);
      [this_out{:}] = feval(fun, s(index));
      if numel(this_out) == 1, this_out = this_out{1}; end
      vout{index} = this_out;
    elseif nout == 0
      feval(fun, s(index));
    else % nout < 0: variable number of output arguments
      this_out = feval(fun, s(index));
      vout{index} = this_out;
    end
  end
  if ~isempty(vout) && numel(vout) == numel(s) && ~isscalar(vout)
    vout = reshape(vout, size(s));
  end

% when all items are of same class and size, make it uniform
  if numel(vout) > 1
    vout_class  = cellfun(@class, vout, 'UniformOutput',false);
    vout_numel  = cellfun(@numel, vout);
    if all(strcmp(vout_class{1}, vout_class)) && all(vout_numel == vout_numel(1))
      % all same class and number of elements
      success = false;
      if size(vout{1},1) == 1 % row vector
        try
          vout = cell2mat(vout);
          success = true;
        end
      end
      if ~success
        try
          vout = cell2mat(vout');
          success = true;
        end
      end
    end
  end
% single returned cell ?
  if iscell(vout) && numel(vout) == 1, vout = vout{1}; end