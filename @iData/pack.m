function b = pack(a)
% b = pack(a) : compress iData storage to save memory
%
%   @iData/pack function to save memory when storing data sets
%   This includes to compare sparse/full storage for large matrices, and
%   reduce the Command history.
%
% input:  s: object or array (iData)
% output: f: compressed object or array (iData)
% ex:     b=pack(a);
%
% See also iData, iData/sparse, iData/full, iData/saveas

if length(a) > 1
  b = a;
  for index=1:length(a(:))
    b(index) = pack(a(index));
  end
  return
end

% extract field type and size
[match, types, nelements]=findfield(a);
b = copyobj(a);

largemat = find(nelements > 1000);
for index=1:length(largemat)
  if ~isempty(strmatch(types{index}, {'double','single','logical'}))
    f = match{index}; % field name
    d = get(a, f);    % content
    if issparse(d),
      who_sparse = whos('d'); 
      d = full(d);
      who_full   = whos('d'); 
    else
      who_full   = whos('d'); 
      d = sparse(d);
      who_sparse = whos('d'); 
    end
    if who_sparse.bytes < who_full.bytes/2, d = sparse(d);
    else                                    d = full(d); end
    set(b, f, d);
  end
end

% now reduce the size of the Command history
h = b.Command;
largemat = cellfun('length', h);
largemat = find(nelements > 1000);
for index=1:length(largemat)
  d = h{index};
  d = [ d(1:800) ' ... ' d(end-100:end) ];
  h{index} = d;
end
b.Command = h;

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),b);
end

