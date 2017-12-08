function [s,lambda,distance,chwidth,energy,wavevector] = Sqw_search_lambda(s)
  % search for the wavelength etc in the object parameters
  
  lambda = []; distance = []; chwidth=[]; energy=[]; wavevector=[];
  if ~isfield(s, 'parameters')
    parameters = [];
  else
    parameters = get(s, 'parameters');
  end
  
  % get lambda, energy, wavevector and check for distance and channel width ----
  % first solution:  search in 'parameters'
  % second solution: findfield in the object
  if isfield(parameters, 'Wavelength'),         
    lambda     = parameters.Wavelength; 
  else
    lambda     = getfieldvalue(s, {'wavelength' 'lambda'});
  end
  if isfield(parameters, 'IncidentEnergy'),     
    energy     = parameters.IncidentEnergy; 
  else
    energy     = getfieldvalue(s, {'IncidentEnergy' 'fixed_energy' 'energy' 'ei'});
  end
  if isfield(parameters, 'IncidentWavevector'), 
    wavevector = parameters.IncidentWavevector; 
  else
    wavevector = getfieldvalue(s, {'IncidentWavevector','wavevector' 'ki'});
  end
  if isfield(parameters,'Distance')
    distance = parameters.Distance;
  else
    distance = getfieldvalue(s, {'Distance_Det_Sample','detector_distance', 'distance'});
  end
  if isfield(parameters, 'ChannelWidth')
    chwidth = parameters.ChannelWidth;
  else
    chwidth = getfieldvalue(s, {'Channel_width','ChannelWidth'});
  end
  
  % now we test if the retyrieved value are OK
  if isempty(lambda)
    if     ~isempty(energy)
      if ischar(energy), energy = get(s, energy); end
      if ~isempty(energy) && energy > 0
        lambda = sqrt(81.805./energy);
      end
    elseif ~isempty(wavevector)
      if ischar(wavevector), wavevector = get(s, wavevector); end
      if ~isempty(wavevector) && wavevector > 0
        lambda = 2*pi./wavevector; 
      end
    end
  end
  
  if isempty(lambda)
    disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' undefined incident neutron wavelength/wavevector/energy.' ]);
    disp('    Define e.g. s.Wavelength=<lambda in Angs> or s.IncidentEnergy=<energy in meV>');
  else
    lambda = mean(lambda(:));
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <wavelength>               =' num2str(lambda) ' [Angs]']);
    setalias(s, 'IncidentWavelength', energy, 'Incident neutron Wavelength [Angs-1]');
  end
  
  if ~isempty(lambda) && lambda > 0 && (isempty(energy) || ischar(energy))
    energy      = 81.805./lambda^2; 
    setalias(s, 'IncidentEnergy', energy, 'Incident neutron Energy [meV]');
  end
  if ~isempty(lambda) && lambda > 0 && (isempty(wavevector)  || ischar(wavevector))
    wavevector  = 2*pi./lambda;
    setalias(s, 'IncidentWavevector', energy, 'Incident neutron Wavevector [Angs-1]');
  end
  
  % search for a sample-to-detector distance
  if ~isempty(distance)
    distance = mean(distance(:));
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <sample-detector distance> =' num2str(distance) ' [m]' ]);
    setalias(s, 'Distance_Sample_Detector', distance);
  end
  
  % search for the Channel Width
  if ~isempty(chwidth)
    chwidth = mean(chwidth);
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <channel width>            =' num2str(chwidth) ' [time unit: s, ms or us]']);
    setalias(s, 'Time_Channel_Width', chwidth);
  end
  
  % ----------------------------------------------------------------------------
  function v = getfieldvalue(s, fields)
    v = [];
    for f=fields(:)'
      if isfield(s, f{1})
        v = get(s, f{1}); return; 
      end
      link = findfield(s, f{1}, 'first exact');
      if isempty(link) && numel(f{1}) > 3
        link = findfield(s, f{1}, 'first');
      end
      if ~isempty(link)
        v = get(s, link); return; 
      end
    end
  
