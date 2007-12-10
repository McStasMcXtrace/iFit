function a=load_ill_tas(a)
% function a=load_ill_tas(a)
%
% Simple postprocessing for ILL/TAS files
%

% get the main data block name
DataBlock        = a.Signal;

% get the main data block header
columns_header   = findstr(a, 'DATA_:','case');
columns_header   = columns_header{1};

% Find spaces and determine proper aliases for the columns
columns = strread(columns_header,'%s','delimiter',' ');

% restrict to the number of columns in DataBlock
c       = size(a, 2);
columns = columns((end-c+1):end);

% compute the normalize variance of each column
index_hkle=[]; % index of QH QK QL EN columns
index_m12 =[]; % index of M1 M2 monitors
Variance = zeros(1,length(columns));
for j=1:length(columns)
  setalias(a,columns{j},a.Signal(:,j));
  if strmatch(columns{j}, {'QH','QK','QL','EN'},'exact')
    index_hkle = [ index_hkle j ];
  end
  if strmatch(columns{j}, {'M1','M2'},'exact')
    index_m12 = [ index_m12 j ];
  end
  if isempty(strmatch(columns{j},{'PNT','CNTS','TI'}, 'exact'))
    if (mean(a.Signal(:,j)))
      Variance(j) = sqrt(sum(a.Signal(:,j).^2)/length(a.Signal(:,j)))/abs(mean(a.Signal(:,j)));
    end
  end
end
% Signal is in CNTS field, 1st axis is probably field with
% remaining greatest variance

% get the monitor
if ~isempty(index_m12)
  [dummy, index]=max(sum(a.Signal(:,index_m12)));
  index_m12 = index_m12(index);
  setalias(a,'Monitor',columns{index_m12},[ 'Monitor ' columns{index_m12} ]);
end
% try with usual scanned variables 'QH','QK','QL','EN'
index = [];
if ~isempty(index_hkle)
  [dummy, index]=max(Variance(index_hkle));
  if dummy > 1e-3
    index = index_hkle(index);
  end
end
if isempty(index)
  [dummy, index]=max(Variance);
end


% retrieve specific information
LOCAL = a.Headers.MetaData.LOCAL; LOCAL=deblank(LOCAL(7:end));
USER  = a.Headers.MetaData.USER;  USER =deblank(USER(7:end));
EXPNO = a.Headers.MetaData.EXPNO; EXPNO=deblank(EXPNO(7:end));
INSTR = a.Headers.MetaData.INSTR; INSTR=deblank(INSTR(7:end));
DATE  = a.Headers.MetaData.DATE;  DATE =deblank(DATE(7:end));
COMND = a.Headers.MetaData.COMND; COMND=deblank(COMND(7:end));
% update object
a.Date = DATE;
a.User = [ EXPNO '=' USER '/' LOCAL '@' INSTR ];
a.Title= [ COMND ' ' a.Title ' ' EXPNO '@' INSTR ];

% make up Signal label

% set Signal and default axis
setalias(a,'Signal','CNTS',[ 'Data CNTS' ]);
setaxis(a,1,columns{index});

