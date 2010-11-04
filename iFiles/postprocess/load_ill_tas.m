function a=load_ill_tas(a)
% function a=load_ill_tas(a)
%
% Simple postprocessing for ILL/TAS files.
% Supports ILL TAS files, including those with multidetectors.
%

a=iData(a);

% re-assign DATA in case the file has more than one data block (MULTI, ...)
try
  DATA  = a.MetaData.DATA; 
  setalias(a,'Signal',DATA,[ 'Data CNTS' ]);
catch
  DATA  = a.Signal;
end

try
  STEPS = a.Data.STEPS;
catch
  STEPS=[];
end

% get the main data block header
[columns_header, data]   = findstr(a, 'DATA_:','case');
if isempty(columns_header), 
  warning(mfilename, [ 'The loaded data set ' a.Tag ' is not an ILL TAS data format' ]);
  return; 
end
columns_header   = columns_header{1};

% detect if this is a MULTI detector file
try
  MULTI  = a.MetaData.MULTI; 
  setalias(a, 'MULTI', 'this.Data.MetaData.MULTI', 'Multidetector');
catch
  MULTI=[];
end
% Find spaces and determine proper aliases for the columns
columns = strread(columns_header,'%s','delimiter',' ;');

% restrict to the number of columns in DataBlock
c       = size(a, 2);
columns = columns((end-c+1):end);

% check if a scan step exists
if ~isempty(STEPS)
  steps_val=struct2cell(STEPS);
  steps_lab=fieldnames(STEPS);
  STEPS=[];
  % get the first non zero step
  for index=1:length(steps_val)
    if abs(steps_val{index}) > 1e-4 && isempty(STEPS)
      STEPS = steps_lab{index};
    end
  end
end
% compute the normalized variance of each column
index_hkle=[]; % index of QH QK QL EN columns
index_m12 =[]; % index of M1 M2 monitors
index_temp=[]; % index for temperatures
index_pal =[]; % index for polarization analysis
Variance = zeros(1,length(columns));
for j=1:length(columns)
  setalias(a,columns{j},a.Signal(:,j)); % create an alias for each column
  % use STEPS
  if (~isempty(strmatch(columns{j}, {'QH','QK','QL','EN'},'exact'))) | ...
    (~isempty(STEPS) & ~isempty(strmatch(columns{j}, STEPS, 'exact'))) | ...
    (~isempty(STEPS) & ~isempty(strmatch([ 'D' columns{j} ], STEPS, 'exact')))
    if isempty(find(index_hkle == j)), index_hkle = [ index_hkle j ]; end
  end
  % and other usual columns
  if strmatch(columns{j}, {'M1','M2'},'exact')
    index_m12 = [ index_m12 j ];
  end
  if strmatch(columns{j}, {'TT','TRT'},'exact')
    index_temp= [ index_temp j ];
  end
  if strmatch(columns{j}, {'PAL'},'exact')
    index_pal= [ index_pal j ];
  end
  if isempty(strmatch(columns{j},{'PNT','CNTS','TI'}, 'exact'))
    if length(a.Signal(:,j))
      Variance(j) = sum( abs(a.Signal(:,j)-mean(a.Signal(:,j)) )) /length(a.Signal(:,j));
    end
  end
end

% Signal is in CNTS field, 1st axis is probably field with
% remaining greatest variance
TT = [];
if ~isempty(index_temp)
  TT = mean(a.Signal(:,index_temp(1)));
else
  index_temp=findfield(a, {'TT','TRT'}, 'case');
  if ~isempty(index_temp), TT = getfield(a, index_temp); end
end

FX = []; KFIX = [];
index_fx  =findfield(a, 'FX', 'case');
if ~isempty(index_fx), FX = getfield(a, index_fx); end
index_kfix=findfield(a, 'KFIX', 'case');
if ~isempty(index_fx), KFIX = getfield(a, index_kfix); end

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
  [dummy, index]=max(Variance); % then set axis as the one that varies the most
end


% retrieve specific information
try
  LOCAL = a.Headers.MetaData.LOCAL; LOCAL=strtrim(LOCAL(7:end));
catch
  LOCAL='';
end
try
  TITLE = a.Headers.MetaData.TITLE; TITLE=strtrim(TITLE(7:end));
catch
  TITLE='';
end
try
  USER  = a.Headers.MetaData.USER;  USER =strtrim(USER(7:end));
catch
  USER='';
end
try
  EXPNO = a.Headers.MetaData.EXPNO; EXPNO=strtrim(EXPNO(7:end));
catch
  EXPNO='';
end
try
  INSTR = a.Headers.MetaData.INSTR; INSTR=strtrim(INSTR(7:end));
catch
  INSTR='';
end
try
  DATE  = a.Headers.MetaData.DATE;  DATE =strtrim(DATE(7:end));
catch
  DATE='';
end
try
  COMND  = a.Headers.MetaData.COMND;  COMND =strtrim(COMND(7:end));
catch
  COMND='';
end

setalias(a, 'COMND', 'this.Data.Headers.MetaData.COMND', 'TAS command');
setalias(a, 'INSTR', 'this.Data.Headers.MetaData.INSTR', 'Instrument used');
setalias(a, 'EXPNO', 'this.Data.Headers.MetaData.EXPNO', 'Experiment number');
setalias(a, 'TITL',  'this.Data.Headers.MetaData.TITLE', 'Dataset title');
% update object
if ~isempty(DATE), a.Date = DATE; end
a.User = [ EXPNO ' ' USER '/' LOCAL '@' INSTR ];
if isempty(TITLE), a.Title= [ COMND ';' a.Title ];
else a.Title= [ TITLE ';' a.Title ]; end
a.Title=regexprep(a.Title,'\s+',' ');


% set Signal and default axis
if isempty(MULTI)
  setalias(a,'Signal','CNTS',[ 'Data CNTS' ]);  % has been defined in DATA_ columns
  setaxis(a,1,columns{index});
else
  setalias(a,'Signal', 'MULTI',[ 'Data CNTS' ]);  % has been defined from MULTI_
  setalias(a,'Channel',1:size(MULTI,2), 'Detector channel');
  setaxis(a,1,columns{index});
  setaxis(a,2,'Channel');
  setalias(a,'Monitor', get(a, columns{index_m12}) * 1:size(MULTI,2));
end
% make up Signal label
xl = xlabel(a);
if ~isempty(TT) & isnumeric(TT), xl = sprintf('%s T=%.2f K',xl,TT); end
if ~isempty(FX) & isnumeric(FX), 
  if ~isempty(KFIX) & isnumeric(KFIX)
    if FX == 1, xl = sprintf('%s Ki=%.2f',xl,KFIX);
    else        xl = sprintf('%s Kf=%.2f',xl,KFIX);
    end
  end
end
xlabel(a, xl);

% handle polarization analysis files
if ~isempty(index_pal)
  % get number of PAL states
  pal = unique(DATA(:,index_pal(1)));
  b = [];
  for j=1:length(pal)
    % create one iData per PAL state
    index_rows = find( DATA(:,index_pal(1)) == pal(j) );
    this_b = a(index_rows, :);
    title(this_b, [ title(this_b) ' [PAL=' num2str(j) ']' ]);
    this_b.Title = [ 'PAL=' num2str(j) ';' this_b.Title  ];
    b = [ b this_b ];
  end
  a = b;
end
  
