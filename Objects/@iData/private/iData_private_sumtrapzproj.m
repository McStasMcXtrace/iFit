function [s, sigma] = iData_private_sumtrapzproj(a,dim, op)
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
% Version: $Date$
% See also iData, iData/plus, iData/prod, iData/cumsum, iData/mean, iData/camproj, iData/trapz

% handle input iData arrays
sigma = [];
if numel(a) > 1
  s = [];
  for index=1:numel(a)
    s = [ s feval(mfilename, a(index), dim, op) ];
  end
  s = reshape(s, size(a));
  return
end

if isscalar(a)
  sigma= double(get(a,'Error'));
  s    = double(a);
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
if isempty(s) || ~isnumeric(s), if dim, s=iData; else s=0; end; return; end
e = get(a,'Error');   % raw Error  (no Monitor weight)
m = iData_private_cleannaninf(get(a,'Monitor'));
if numel(e) > 1 && all(e(:) == e(1)), e=e(1); end
if numel(m) > 1 && all(m(:) == m(1)), m=m(1); end
if isempty(m), m=1; end

% take into account the Monitor
if not(all(m(:) == 0 | m(:) == 1))
  s = genop(@rdivide, s, m); e = genop(@rdivide,e,m); 
end

s = iData_private_cleannaninf(s);
e = iData_private_cleannaninf(e);  

[link, lab] = getalias(a, 'Signal');
cmd= a.Command;
b  = copyobj(a);

if any(dim == 0) && ~strcmp(op, 'sum')
  dim=1:ndims(a);
end

myisvector = @(c)length(c) == numel(c);

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
      if dim(index) == 1 && myisvector(s)
        s = s(:); e=e(:); m=m(:);
      end
      if numel(e) > 1, e = feval(op, e, dim(index)); end % sum(e2/m2)
      if numel(m) > 1, m = feval(op, m, dim(index)); end % sum(m)
      s = feval(op, s, dim(index));                      % sum(s/m)
    end
    % Store Signal
    s=squeeze(s); e=sqrt(squeeze(e)); m=squeeze(m);
    if not(all(m(:) == 0 | m(:) == 1)), s=genop(@times,s,m); e=genop(@times,e,m); end
    setalias(b,'Signal', s, [op ' of ' lab ' along axis ' num2str(dim) ]);
    
  case {'prod','cumprod'} % PROD ===============================================
    % product on all requested dimensions 
    for index=1:numel(dim)
      if dim(index) == 1 && myisvector(s)
        s = s(:); e=e(:); m=m(:);
      end
      if numel(e) > 1, e = feval(op, s+e/2, dim(index))-feval(op, s-e/2, dim(index)); end
      if numel(m) > 1, m = feval(op, m, dim(index)); end
      s = feval(op, s, dim(index)); 
    end
    % Store Signal
    s=squeeze(s); e=squeeze(e); m=squeeze(m);
    if not(all(m(:) == 0 | m(:) == 1)), s=genop(@times,s,m); e=genop(@times,e,m); end
    setalias(b,'Signal', s, [op ' of ' lab ' along axis ' num2str(dim) ]);
    
  case {'trapz','cumtrapz'} % TRAPZ ============================================

    % integral is: trapz(Signal/Monitor)              =trapz(s)
    % error bar:   sigma=sqrt(trapz(Error2/Monitor2)) =sqrt(trapz(e2))
    
    e=e.^2;

    for index=1:numel(dim)
      [x, xlab]     = getaxis(a,dim(index));

      if dim(index) ~= 1  % we put the dimension to integrate on as first
        perm=1:ndims(a);
        perm(dim(index))=1; perm(1)=dim(index);
        s = permute(s, perm);
        e = permute(e, perm); 
        m = permute(m, perm);
      elseif myisvector(s)
        s = s(:); x=x(:); e=e(:); m=m(:);
      end
      % make the integration
      if ~isscalar(s)
        % check dimension of axis x
        if numel(x) == prod(size(x)), x=x(:); end
        if numel(x) ~= size(s, 1), x=x(:,1); end
        if numel(e) > 1, e = feval(op, x, e, 1); end % trapz(x,e2/m2,1)
        if numel(m) > 1, m = feval(op, 1:length(x), m, 1); end % sum(m)
        % check length(x) == size(s,1)
        if numel(x) ~= size(s,1)
           x0 = unique(x);
           if numel(x0) == size(s,1)
               x = x0;
               a = setaxis(a, dim(index), x);
           end
        end
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
    setalias(b,'Signal', s, [ op ' of ' lab ' along ' xlab ]);

  case 'camproj' % camproj =====================================================
    % accumulates on all axes except the rank specified
    [x, xlab]     = getaxis(a,dim);
    s             = subsref(a,struct('type','.','subs','Signal')); 
    if myisvector(s),   lx=length(s); else 
      lx=size(s, dim); 
      if myisvector(x)
        % replicate X along other axes
        sz = size(s); sz(dim) = 1;
        x = repmat(x, sz);
      end
    end
    rmaxis(b);
    setaxis(b, 1, x(:));
    % Store Signal
    setalias(b,'Signal', s(:), [ 'projection of ' lab ' on axis ' num2str(dim) ]);     % Store Signal
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
        % x = getaxis(a, index); 
        [x, xlab] = getaxis(a, num2str(index)); % get axis definition and label
        if ~isnumeric(x), x = getaxis(a, index); end
        S.type = '()';
        S.subs = {};
        if ndims(b) == ndims(x)
          S.subs={ ':' };
        else
          for j=1:ndims(a), 
            if j ~= index, S.subs{j}= 1;
            else           S.subs{j}=':'; end
          end
        end
        if ~myisvector(x),  x = subsref(x,S); end
        setaxis(b, ax_index, x);
        label(  b, ax_index, xlab);
        ax_index = ax_index+1;
      end
    end
  end
	
elseif dim == 0
  s = sum(s(:));
  sigma = double(sqrt(sum(e(:).^2)));
  s = double(s);
  return  % scalar
end
b.Command=cmd;
b = iData_private_history(b, op, b, dim);
s = b;

if isscalar(s)
  sigma = double(get(s,'Error'));
  s=double(s); 
end

% reset warnings
iData_private_warning('exit',mfilename);

