function b=load_ill_inx(a)
% function a=load_ill_inx(a)
%
% Returns an iData style dataset from an ILL INX file
%

% handle input iData arrays
if length(a(:)) > 1
  b = [];
  for index=1:length(a(:))
    b(index) = feval(mfilename, a(index));
  end
  return
end

% the INX data set is loaded in catenate mode. 
% As there is not header, the first block is Data_0
% its number of rows is the block number
% its last column is the length of each group

% get list of Data fields and their size
nblocks = size(a.Data.Data_0, 1);
b = [];
blockstart=1;
alldata=struct2cell(a.Data);
blocksdouble=cellfun('isclass', alldata, 'double');
alldata=alldata(find(blocksdouble));
blockslength=alldata{1};
blocksangle =alldata{end-1};
blocksdata  =alldata{end};
for index=1:nblocks
  blocklength=blockslength(index, end);  
  data = [];
  data.group = index;
  data.block = blocksdata((blockstart+1):(blockstart+blocklength), :);  % skip first row ;
  data.angle = blocksangle(index, 1);
  data.wavelength=  blocksangle(index, 2);
  data.wavevector=  blocksangle(index, 3);
  data.temperature= blocksangle(index, 4);
  data.Headers=[];
  c=a; setalias(c, getalias(c)); % copy a and clear Aliases
  c.Data=data;
  c=iData(c);
  c.Title   =[ '#' num2str(index) ' angle=' num2str(data.angle) ' T=' num2str(data.temperature) ' lambda=' num2str(data.wavelength) ' ' a.Title ];
  
  setalias(c,'Signal', 'Data.block(:,2)',[ 'Signal #' num2str(index) ' angle=' num2str(data.angle)  ]);
  setalias(c,'Error',  'Data.block(:,3)');
  setalias(c,'Energy', 'Data.block(:,1)',[ 'Energy [meV] T=' num2str(data.temperature) ' lambda=' num2str(data.wavelength) ]);
  setalias(c,'Angle', 'Data.angle','Detection angle [deg]');
  setalias(c,'Temperature', 'Data.temperature','Sample Temperature [K]');
  setalias(c,'Wavelength', 'Data.wavelength','Incident wavelength [Angs]');
  setaxis(c, 1, 'Energy');
  setaxis(c, 2, 'Angle');
  
  blockstart = blockstart+1+blocklength;
  b = [ b c ];
end



