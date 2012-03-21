function b=load_ill_inx(a)
% function a=load_ill_inx(a)
%
% Returns an iData style dataset from an ILL INX file
%

if ~isa(a,'iData')
  b = load(iData,a,'ILL INX');
  return
end

% handle input iData arrays
if length(a(:)) > 1
  b = [];
  for index=1:length(a(:))
    b(index) = feval(mfilename, a(index));
  end
  return
end

% the data read with read_inx comes as:
% s.Data.header: char
% s.Data.Par:    double with parameters after the header/comment line
% s.Data.Mat:    double
Data = a.Data;
Data.angle       = a.Data.Par(:,1);
Data.wavelength  = a.Data.Par(:,2);
Data.wavevector  = a.Data.Par(:,3);
Data.temperature = a.Data.Par(:,4);
Data.signal      = squeeze(a.Data.Mat(:,2,:));
Data.error       = squeeze(a.Data.Mat(:,3,:));
Data.energy      = squeeze(a.Data.Mat(:,1,:));
Data.energy      = Data.energy(:,1);
a.Data = Data;

b = a;

setalias(b,'Signal', 'Data.signal', b.Data.header(2,:,1));
setalias(b,'Error',  'Data.error');
setalias(b,'Energy', 'Data.energy', ...
  [ 'Energy [meV] T=' num2str(mean(b.Data.temperature)) ' lambda=' num2str(mean(b.Data.wavelength)) ]);
setalias(b,'Angle',       'Data.angle','Angle [deg]');
setalias(b,'Wavelength',  'Data.wavelength','Wavelength [Angs]');
setalias(b,'Wavevector',  'Data.wavevector','Wavevector [Angs]');
setalias(b,'Temperature', 'Data.temperature','Sample Temperature [K]');

if ndims(b) == 1
  setaxis(b, 1, 'Energy');
elseif ndims(b) == 2
  setaxis(b, 1, 'Energy');
  setaxis(b, 2, 'Angle');
end
b.Title   =[ b.Data.header(2,:,1) ' ' a.Title ];

