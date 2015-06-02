function s = iData_private_sumtrapzproj(a,dim, op)
% s = iData_private_sumtrapzproj(a,dim, op) : computes the sum/trapz/camproj of iData objects elements
%
%   @iData/iData_private_sumtrapzproj function to compute the sum/trapz/camproj of the elements of the data set
%
% input:  a: object or array (iData/array of)
%         dim: dimension to accumulate (int/array of/'radial')
%         op: 'sum','trapz','camproj','prod'
% output: s: sum/trapz/camproj of elements (iData/scalar)
% ex:     c=iData_private_sumtrapzproj(a, dim, 'sum');
%
% Version: $Revision$
% See also iData, iData/plus, iData/prod, iData/cumsum, iData/mean, iData/camproj, iData/trapz

% handle input iData arrays
if numel(a) > 1
  s = [];
  for index=1:numel(a)
    s = [ s feval(op, a(index), dim) ];
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

% some cases:
% signal is vector nD -> 

% in all cases except projection, resample the data set on a grid
% make axes single vectors for sum/trapz/... to work
if all(dim > 0) & any(strcmp(op, {'sum','cumsum','prod','cumprod','trapz','cumtrapz'}))
  a = meshgrid(a, 'vector');
end

s = get(a,'Signal');  % raw Signal (no Monitor weight)
e = get(a,'Error');   % raw Error  (no Monitor weight)
m = iData_private_cleannaninf(get(a,'Monitor'));

% take into account the Monitor
if not(all(m(:) == 0 | m(:) == 1))
  s = genop(@rdivide, s, m); e = genop(@rdivide,e,m); 
end

s = iData_private_cleannaninf(s);
e = iData_private_cleannaninf(e);  

[link, label] = getalias(a, 'Signal');
cmd= a.Command;
b  = copyobj(a);

if any(dim == 0) && ~strcmp(op, 'sum')
  dim=1:ndims(a);
end

% in the following, we have e.g. e=e/m; s=s/m;
if all(dim > 0)
  % compute new object: op(Signal/Monitor) op(Error/Monitor)
  switch op
  case {'sum','cumsum'} % SUM ==================================================
    % sum on all requesteddimensions 

    % integral is: sum(Signal/Monitor)              =sum(s)
    % error bar:   sigma=sqrt(sum(Error2/Monitor2)) =sqrt(sum(e2))
    e=e.^2; 
    for index=1:numel(dim)
      if dim(index) == 1 && isvector(s)
        s = s(:); e=e(:); m=m(:);
      end
      if numel(e) > 1, e = feval(op, e, dim(index)); end % sum(e2/m2)
      if numel(m) > 1, m = feval(op, m, dim(index)); end % sum(m)
      s = feval(op, s, dim(index));                      % sum(s/m)
    end
    % Store Signal
    s=squeeze(s); e=sqrt(squeeze(e)); m=squeeze(m);
    if not(all(m(:) == 0 | m(:) == 1)), s=genop(@times,s,m); e=genop(@times,e,m); end
    setalias(b,'Signal', s, [op ' of ' label ' along axis ' num2str(dim) ]);
    
  case {'prod','cumprod'} % PROD ===============================================
    % product on all requested dimensions 
    for index=1:numel(dim)
      if dim(index) == 1 && isvector(s)
        s = s(:); e=e(:); m=m(:);
      end
      if numel(e) > 1, e = feval(op, s+e/2, dim(index))-feval(op, s-e/2, dim(index)); end
      if numel(m) > 1, m = feval(op, m, dim(index)); end
      s = feval(op, s, dim(index)); 
    end
    % Store Signal
    s=squeeze(s); e=squeeze(e); m=squeeze(m);
    if not(all(m(:) == 0 | m(:) == 1)), s=genop(@times,s,m); e=genop(@times,e,m); end
    setalias(b,'Signal', s, [op ' of ' label ' along axis ' num2str(dim) ]);
    
  case {'trapz','cumtrapz'} % TRAPZ ============================================

    % integral is: trapz(Signal/Monitor)              =trapz(s)
    % error bar:   sigma=sqrt(trapz(Error2/Monitor2)) =sqrt(trapz(e2))
    
    e=e.^2;

    for index=1:numel(dim)
      [x, xlab]     = getaxis(a,dim(index)); x=x(:);

      if dim(index) ~= 1  % we put the dimension to integrate on as first
        perm=1:ndims(a);
        perm(dim(index))=1; perm(1)=dim(index);
        s = permute(s, perm);
        e = permute(e, perm); 
        m = permute(m, perm);
      elseif isvector(s)
        s = s(:); x=x(:); e=e(:); m=m(:);
      end
      % make the integration
      if ~isscalar(s)
        if numel(e) > 1, e = feval(op, x, e, 1); end % trapz(x,e2/m2,1)
        if numel(m) > 1, m = feval(op, 1:length(x), m, 1); end % sum(m)
        s = feval(op, x, double(s), 1); % trapz(x,s,1)

        if dim(index) ~= 1  % restore initial axes
          s = permute(s,perm);
          e = permute(e,perm);
          m = permute(m,perm);
        end
      end
    end

    % Store Signal
    s=squeeze(s); e=sqrt(squeeze(e)); m=squeeze(m);
    if not(all(m(:) == 0 | m(:) == 1)), s=genop(@times,s,m); e=genop(@times,e,m); end
    setalias(b,'Signal', s, [ op ' of ' label ' along ' xlab ]);

  case 'camproj' % camproj =====================================================
    % accumulates on all axes except the rank specified
    [x, xlab]     = getaxis(a,dim);
    s             = subsref(a,struct('type','.','subs','Signal')); 
    if isvector(s),   lx=length(s); else 
      lx=size(s, dim); 
      if isvector(x)
        % replicate X along other axes
        sz = size(s); sz(dim) = 1;
        x = repmat(x, sz);
      end
    end
    rmaxis(b);
    setaxis(b, 1, x(:));
    % Store Signal
    setalias(b,'Signal', s(:), [ 'projection of ' label ' on axis ' num2str(dim) ]);     % Store Signal
    b = setalias(b, 'Error', abs(e(:)));
    b = setalias(b, 'Monitor', m(:));
    b = hist(b, lx); % faster than doing sum on each dimension

  end % switch (op, compute)
  
  % store new object
  b = setalias(b, 'Error', abs(e));
  b = setalias(b, 'Monitor', m);
	
  if any(strcmp(op, {'sum','trapz','prod'}))
    rmaxis(b);  % remove all axes, will be rebuilt after operation
    % put back initial axes, except those integrated
    ax_index=1;
    for index=1:ndims(a)
      if all(dim ~= index)
        [x, xlab] = getaxis(a, num2str(index)); % get axis definition and labelget(a
        setaxis(b, ax_index, x);
        ax_index = ax_index+1;
      end
    end
  end
	
elseif dim == 0
  s = sum(s(:));
  s = double(s);
  return  % scalar
end
b.Command=cmd;
b = iData_private_history(b, op, b, dim);
s = b;

if isscalar(s), s=double(s); end

% reset warnings
iData_private_warning('exit',mfilename);

