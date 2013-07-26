function [Data, Att] = read_hdf5(hdf_fname)
% HDF5 Extract read and return HDF5 data.
% DATASTRUCT = HDF5EXTRACT('FILENAME')

% by Daniel Buckton
% 28 Mar 2007 (Updated 28 Mar 2007) 

DataInfo = hdf5info(hdf_fname);
BaseStr = 'DataInfo.GroupHierarchy';

LevelStr.CurrentLevel =1;
LevelStr.MaxLevel=16;
Data = []; Att = [];
[Data, Att] = get_struct(Data,Att, DataInfo,hdf_fname,BaseStr,LevelStr);
Data.Attributes = Att;

function [Data, Att] = get_struct(Data,Att,DataInfo,hdf_fname,DataStr,LevelStr)

if LevelStr.CurrentLevel > LevelStr.MaxLevel
    error('Exceeded Maximum Level')    
end

nGroups     = eval(['length(' DataStr '.Groups)']);
nDataSets   = eval(['length(' DataStr '.Datasets)']);
nAttributes = eval(['length(' DataStr '.Attributes)']);

for iDataSet = 1:nDataSets;
  % get Data and Attributes
    eval( [ '[D,A] = hdf5read(''' hdf_fname ''',''' eval([DataStr '.Datasets(' num2str(iDataSet) ').Name']) ''', ''ReadAttributes'',true);' ] )
    if strcmp(class(D), 'hdf5.h5string'), D=char(D.Data); end
    eval(['Data.' eval([ 'ProcessString(' DataStr '.Datasets(' num2str(iDataSet) ').Name)']) ' = D;' ]); clear D;
    if ~isempty(A)
      eval(['Att.' eval([ 'ProcessString(' DataStr '.Datasets(' num2str(iDataSet) ').Name)']) ' = A;' ])
    end
end

LevelStr.CurrentLevel =LevelStr.CurrentLevel +1;
for iGroup = 1:nGroups,
    DataInputStr=[DataStr '.Groups(' num2str(iGroup) ')'];
    [Data,Att] = get_struct(Data,Att,DataInfo,hdf_fname,DataInputStr,LevelStr);
end


function s1 = ProcessString(s)
%PROCESSSTRING Replace matlab incompatable strings with underscore or .
%when a tree structure is desired
%   PROCESSSTRING(S) replaces incompatable characters in string S with _.
%

if isempty(s)
   s1 = s([]);
else
   
   if ~isstr(s),
      warning('MATLAB:deblank:NonStringInput','Input must be a string.')
   end
   
   s1 = s;
   
   % Replace incompatable characters
   ind = find( isspace(s) | (s =='-') | (s =='=') | (s =='.'));
   if ~isempty(ind),
      if ind(1)==1,
          s1(1) = [];
          ind(1) = [];
          ind = ind-1;
      end
      
      s1(ind) = repmat('_',1,length(ind));
   end

   s = s1;
   
   ind = find( s =='/');
   if ~isempty(ind),
      if ind(1)==1,
          s1(1) = [];
          ind(1) = [];
          ind = ind-1;
      end
      s1(ind) = repmat('.',1,length(ind));
      ind2 = find(s1(ind+1)=='_');
      s1(ind(ind2)+1)=[];
  end
end


function s1 = defalse(s)
%DEFALSE Replace matlab incompatable strings with underscore.
%   DEFALSE(S) replaces incompatable characters in string S with _.
%

if isempty(s)
   s1 = s([]);
else
   
   if ~isstr(s),
      warning('MATLAB:deblank:NonStringInput','Input must be a string.')
   end
   
   s1 = s;
   
   % Replace incompatable characters
   ind = find( isspace(s) | (s == '/') | (s =='-') | (s =='=') | (s =='.'));
   if ~isempty(ind),
      if ind(1)==1,
          s1(1) = [];
          ind(1) = [];
          ind = ind-1;
      end
      
      s1(ind) = repmat('_',1,length(ind));
   end
   if s1(1)=='_'
       s1(1)=[];
   end
end

function s1 = endstring(s)
%ENDSTRING Returns any part of the string that is after a foreslash (/).

if isempty(s)
   s1 = s([]);
else
   
   if ~isstr(s),
      warning('MATLAB:endstring:NonStringInput','Input must be a string.')
   end
   
   s1 = s;
   
   if sum(s == '/')~=0,
      ind = find(s(end:-1:1)=='/');
      s1 = s1(end-ind(1)+2:end);       
   end
             
end
