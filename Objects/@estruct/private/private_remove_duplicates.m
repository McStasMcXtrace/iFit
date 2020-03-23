function out = private_remove_duplicates(out)
% PRIVATE_REMOVE_DUPLICATES Check for duplicate data sets.

% clean 'out' for unique entries (e.g. mcstas/mcstas.sim creates duplicates)
[b, i1] = unique(get(out, 'Source')); % name is unique
if 1 < length(out) && length(i1) < length(out)
  % some data sets seem to be duplicated: make additional tests
  % look for similarities
  sources = get(out, 'Source');
  titls   = get(out, 'Title');
  labs    = get(out, 'Label');
  sums    = sum(out, 0); % total signal
  i       = 1:length(out);
  for index=1:length(out)
    j = find(strcmp(sources{index}, sources(:)) & strcmp(titls{index}, titls(:)) & strcmp(labs{index}, labs(:)) & sums(index) == sums(:));
    if length(j) > 1, i(j(2:end)) = 0; end
  end
  removed = find(i == 0);
  i       = unique(i(i>0));
  if length(out) > length(i)
    if out(1).verbose
      warning('%s: Removing duplicated data sets %i -> %i', mfilename, length(out), length(i))
      warning([ mfilename ': ' char(out(removed)) ])
    end
    out = out(i);
  end
end
