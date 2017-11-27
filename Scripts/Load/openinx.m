function out = openinx(filename)
%OPENINX Open an INX tof data file, display it
%        and set the 'ans' variable to an iData object with its content
% (c) E.Farhi, ILL. License: EUPL.

if ~isa(filename,'iData')
  out = iData(iLoad(filename,'ILL INX')); % no post-processing
else
  out = filename;
end
clear filename;

if numel(out) > 1
  % handle input iData arrays
  for index=1:numel(out)
    out(index) = feval(mfilename, out(index));
  end
elseif isfield(out.Data, 'header') && isfield(out.Data, 'Par') && isfield(out.Data, 'Mat')
  % the data read with read_inx comes as:
  % s.Data.header: char
  % s.Data.Par:    double with parameters after the header/comment line
  %   Angle, incident energy, transfered wave-vector (Q), mean temperature, -, -, -, channel width (musec), -
  Data = out.Data;
  Data.angle              = Data.Par(:,1);
  Data.IncidentEnergy     = Data.Par(:,2);
  Data.IncidentWavevector = Data.Par(:,3);
  Data.temperature        = Data.Par(:,4);
  Data.signal      = squeeze(Data.Mat(:,2,:));
  Data.error       = squeeze(Data.Mat(:,3,:));
  Data.wavelength  = sqrt(81.805./Data.IncidentEnergy);
  % a few checks for consistency
  if abs(2*pi/mean(Data.wavelength) - Data.IncidentWavevector) > 0.01*Data.IncidentWavevector
    disp([ mfilename ': WARNING: The loaded data set ' out.Tag ' from ' out.Source ' has an inconsistent Energy:' ])
    disp(sprintf('  Incident Energy=%g [meV] Wavevector=%g [Angs-1] Wavelength=%g [Angs]', ...
      mean(Data.IncidentEnergy), mean(Data.IncidentWavevector), mean(Data.wavelength)));
  end
  if ~Data.angle  % this is an Angle data set, not Energy
    isAngleData = 1;
    Data.angle     = squeeze(Data.Mat(:,1,:));
    Data.angle     = Data.angle(:,1);
    Data.energy    = 81.805/Data.wavelength^2;
  else
    isAngleData = 0;
    Data.energy      = squeeze(Data.Mat(:,1,:));
    Data.energy      = Data.energy(:,1);
  end
  if numel(unique(Data.temperature)) == 1, Data.temperature=Data.temperature(1); end
  if numel(unique(Data.wavelength))  == 1, Data.wavelength =Data.wavelength(1); end
  if numel(unique(Data.IncidentWavevector))  == 1, Data.IncidentWavevector =Data.IncidentWavevector(1); end
  if numel(unique(Data.IncidentEnergy))  == 1,     Data.IncidentEnergy     =Data.IncidentEnergy(1); end
  out.Data = Data; 

  setalias(out,'Signal', 'Data.signal', out.Data.header(2,:,1));
  setalias(out,'Error',  'Data.error');
  if isAngleData
    setalias(out,'Angle', 'Data.angle', ...
      [ 'Angle [deg] T=' num2str(mean(out.Data.temperature)) ' lambda=' num2str(mean(out.Data.wavelength)) ]);
  else
    setalias(out,'Energy', 'Data.energy', ...
      [ 'Energy [meV] T=' num2str(mean(out.Data.temperature)) ' lambda=' num2str(mean(out.Data.wavelength)) ]);
    setalias(out,'Angle',       'Data.angle','Angle [deg]');
  end
  setalias(out,'IncidentEnergy',      Data.IncidentEnergy,    'Incident Energy [meV]');
  setalias(out,'Wavelength',          Data.wavelength,        'Incident Wavelength [Angs]');
  setalias(out,'IncidentWavevector',  Data.IncidentWavevector,'Incident Wavevector [Angs-1]');
  setalias(out,'Temperature', 'Data.temperature','Sample Temperature [K]');

  if ndims(out) == 1
    if isAngleData
      setaxis(out, 1, 'Angle');
    else
      setaxis(out, 1, 'Energy');
    end
  elseif ndims(out) == 2
    setaxis(out, 1, 'Energy');
    setaxis(out, 2, 'Angle');
  end
  out.Title   =[ out.Data.header(2,:,1) ' ' out.Title ];
else
  warning([ mfilename ': The loaded data set ' out.Tag ' from ' out.Source ' is not an ILL INX data format.' ]);
end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end
