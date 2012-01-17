function a=load_psi_RITA(a0)
% function a=load_fig(a0)
%
% Returns an iData style dataset from a PSI RITA2 file (HDF5)
%
% Version: $Revision: 1.1 $
% See also: iData/load, iLoad, save, iData/saveas

% handle input iData arrays
if length(a0(:)) > 1
  for index=1:length(a0(:))
    a(index) = feval(mfilename, a0(index));
  end
  return
end

if isempty(findfield(a0, 'RITA_2'))
  a = a0;
  return
end

% RITA2 file identified...
a = a0;
% Alias a.Data.entry1.data members
data_alias = fieldnames(a.Data.entry1.data);
for index=1:length(data_alias)
  setalias(a, data_alias{index}, [ 'Data.entry1.data.' data_alias{index} ]);
end

% Alias a.Data.entry1.RITA_2.detectorwindows.counts
setalias(a, 'windowcounts', 'Data.entry1.RITA_2.detectorwindows.counts');
setalias(a, 'Signal',       'Data.entry1.data.counts');

% now we create 9 objects per window: 
% a.Data.entry1.RITA_2.detectorwindows.counts(1:9, :)
% a.Data.entry1.data -> 2 dimensional members will extract (1:9,:)
b = zeros(a, [ 1 9 ]);
index_axis    = 1;
for j=1:9
  this = b(j);
  index_axis_std= 0;
  for index=1:length(data_alias)
    member=get(this, [ 'Data.entry1.data.' data_alias{index} ]);
    if ndims(member)== 2 && ~isvector(member)
      setalias(this, data_alias{index}, ...
        [ 'Data.entry1.data.' data_alias{index} '(' num2str(j), ',:)' ]);
    end
    % find the member which varies most
    if ~isvector(member), member=member(j,:); end
    member=double(member);
    if ndims(member)==2 && std(member) > index_axis_std && j==5
      if any(strcmp(data_alias{index}, {'Qh','Qk','Ql','energy'}))
        index_axis_std = std(member);
        index_axis = index;
      elseif any(strcmp(data_alias{index}, {'a1','a2','a3','a4','a5','a6'})) && index_axis_std==0
        index_axis_std = std(member);
        index_axis = index;
      end
    end
  end
  setalias(this, 'windowcounts', ...
    [ 'Data.entry1.RITA_2.detectorwindows.counts(' num2str(j) ',:)' ]);
  setalias(this, 'Signal', 'windowcounts');
  this.Title = [ '#' num2str(j) ' ' this.Title ];
  setaxis(this, 1, data_alias{index_axis});
  
  b(j) = this;
end

a = [ a b ];

