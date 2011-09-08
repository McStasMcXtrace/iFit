function a=load_vitess_2d(a)
% function a=load_vitess_2d(a)
%
% Returns an iData style dataset from a VITESS 2d monitor file
%
a=iData(a);

% Vitess 2D have:
% * 2 numerical blocks
% * a first block with one column less that second block
% * first block has one or 2 rows
[match, types, nelements] = findfield(a);
index = strmatch('double',types);
if length(index) ~= 2, return; end

f1 = get(a,match{index(1)});
f2 = get(a,match{index(2)});

if size(f1,2) ~= size(f2,2) -1, return; end
if size(f1,1) > 2,              return; end

x  = f1(1,:);
y  = f2(:,1);
f2 = f2(:, 2:end);

a{1}=y;
a{2}=x;
a{0}=f2;
