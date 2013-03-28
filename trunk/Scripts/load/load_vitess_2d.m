function a=load_vitess_2d(a)
% function a=load_vitess_2d(a)
%
% Returns an iData style dataset from a VITESS 2d monitor file
%
% Version: $Revision$
% See also: iData/load, iLoad, save, iData/saveas

if ~isa(a,'iData')
  a = load(iData,a,mfilename);
  return
end

% handle input iData arrays
if length(a(:)) > 1
  for index=1:length(a(:))
    a(index) = feval(mfilename, a(index));
  end
  return
end

a=iData(a);

% Vitess 2D have:
% * at least 2 numerical blocks, and we use the last two
% * a first block with one column less that second block
% * first block has one or 2 rows
[match, types, nelements] = findfield(a);
index = strmatch('double',types);
if length(index) < 2, return; end

f1 = get(a,match{index(end-1)});
f2 = get(a,match{index(end)});

if size(f1,2) ~= size(f2,2) -1, return; end % 2nd block has more columns that 1st
if size(f1,1) > 2,              return; end % more than 2 lines

a{0}=[ match{index(end)} '(:,2:end)' ];
a{1}=[ match{index(end)} '(:,1)' ];
a{2}=[ match{index(end-1)} '(1,:)' ];
