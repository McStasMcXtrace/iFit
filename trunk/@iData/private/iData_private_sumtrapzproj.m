function s = iData_private_sumtrapzproj(a,dim, op)
% s = iData_private_sumtrapzproj(a,dim, op) : computes the sum/trapz/camproj of iData objects elements
%
%   @iData/iData_private_sumtrapzproj function to compute the sum/trapz/camproj of the elements of the data set
%
% input:  a: object or array (iData/array of)
%         dim: dimension to accumulate (int/array of)
%         op: 'sum','trapz','camproj','prod'
% output: s: sum/trapz/camproj of elements (iData/scalar)
% ex:     c=iData_private_sumtrapzproj(a, dim, 'sum');
%
% Version: $Revision: 1.8 $
% See also iData, iData/plus, iData/prod, iData/cumsum, iData/mean, iData/camproj, iData/trapz

% handle input iData arrays
if numel(a) > 1
  if dim ~= 0, s = a; else s=zeros(size(a)); end
  for index=1:numel(a)
    s(index) = feval(op, a(index), dim);
  end
  s = reshape(s, size(a));
  return
end

if isscalar(a)
  s = double(a);
  return
end

% removes warnings
iData_private_warning('enter',mfilename);

% in all cases, resample the data set on a grid
%a = interp(a,'grid');

S.type = '()';
S.subs = cell(1, ndims(a));
rebin  = 0;
% make axes single vectors for sum/trapz/... to work
for index=1:ndims(a)
  [x, xlab] = getaxis(a, index);
  [dummy, i,j] = unique(x);
  if length(x) < length(i)
    S.subs{index} = i;
    rebin = 1;
  else
    S.subs = ':';
  end
end
if rebin
  a = subsref(a, S);
end

s = iData_private_cleannaninf(get(a,'Signal'));
e = iData_private_cleannaninf(get(a,'Error'));
m = iData_private_cleannaninf(get(a,'Monitor'));

[link, label] = getalias(a, 'Signal');
cmd= a.Command;
b  = copyobj(a);
rmaxis(b);  % remove all axes, will be rebuilt after operation


if all(dim > 0)
  % compute new object
  switch op
  case 'sum' % SUM =============================================================
    % sum on all dimensions requested
    for index=1:length(dim(:))
      if numel(e) > 1, e = sum(e, dim(index)); end
      if numel(m) > 1, m = sum(m, dim(index)); end
      s = sum(s, dim(index)); 
    end
    % Store Signal
    s=squeeze(s); e=squeeze(e); m=squeeze(m);
    setalias(b,'Signal', s, [op ' of ' label ' along axis ' num2str(dim) ]);
    
  case 'prod' % PROD ===========================================================
    % product on all dimensions requested
    for index=1:length(dim(:))
      if numel(e) > 1, e = prod(s+e/2, dim(index))-prod(s-e/2, dim(index)); end
      if numel(m) > 1, m = prod(m, dim(index)); end
      s = prod(s, dim(index)); 
    end
    % Store Signal
    s=squeeze(s); e=squeeze(e); m=squeeze(m);
    setalias(b,'Signal', s, [op ' of ' label ' along axis ' num2str(dim) ]);
    
  case 'trapz' % TRAPZ =========================================================
    for index=1:length(dim(:))
      [x, xlab]     = getaxis(a,dim(index));
      if dim(index) ~= 1  % we put the dimension to integrate on as first
        perm=1:ndims(a);
        perm(dim(index))=1; perm(1)=dim(index);
        s = permute(s, perm);
        e = permute(e, perm); 
        m = permute(m, perm);
      end
      % make the integration
      if numel(e) > 1, e = trapz(x, e); end
      if numel(m) > 1, m = trapz(x, m); end
      s = trapz(x, s);
      if dim(index) ~= 1  % restore initial axes
        s = permute(s,perm);
        e = permute(e,perm);
        m = permute(m,perm);
      end
    end
    % Store Signal
    s=squeeze(s); e=squeeze(e); m=squeeze(m);
    setalias(b,'Signal', s, [ op ' of ' label ' along ' xlab ]); 
    
  case 'camproj' % camproj =====================================================
    % accumulates on all axes except the rank specified
    for index=1:ndims(a)
      if index~=dim, 
        if numel(e) > 1, e = sum(e, index); end
        if numel(m) > 1, m = sum(m, index); end
        s = trapz(s, index); 
      end
    end
    % Store Signal
    setalias(b,'Signal', s, [ 'projection of ' label ' on axis ' num2str(dim) ]);     % Store Signal
	  
	end % switch (compute)
	
	% store new object
	b = set(b, 'Error', abs(e), 'Monitor', m);
	switch op
	case {'sum','trapz','prod'}
    % put back initial axes, except those integrated
    ax_index=1;
    for index=1:ndims(a)
      if all(dim ~= index)
        [x, xlab] = getaxis(a, num2str(index)); % get axis definition and label
        setaxis(b, ax_index, x);
        ax_index = ax_index+1;
      end
    end
  case 'camproj'
    % set projection axis
    [x, xlab] = getaxis(a, num2str(dim)); % get axis definition and label
    setaxis(b, 1, x);

    if dim == 1, 
	    s=size(b);
	    if s(1) == 1, b=transpose(b); end
	  end
	end % switch (store)
	
elseif dim == 0
  if strcmp(op, 'camproj')
    s = feval('sum', a, 1:ndims(a));
  else
    s = feval(op, a, 1:ndims(a));
  end
  s = double(s);
  return  % scalar
end
b.Command=cmd;
b = iData_private_history(b, op, b, dim);
s = b;

% reset warnings
iData_private_warning('exit',mfilename);

