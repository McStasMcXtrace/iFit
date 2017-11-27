function [s,lambda,distance,chwidth,energy,wavevector] = Sqw_search_lambda(s)
  % search for the wavelength etc in the object parameters
  
  lambda = []; distance = []; chwidth=[]; energy=[]; wavevector=[];
  if ~isfield(s, 'parameters'), return; end
  parameters = get(s, 'parameters');
  
  % get lambda, energy, wavevector and check for distance and channel width
  if isfield(parameters, 'Wavelength'),         lambda     = parameters.Wavelength; end
  if isfield(parameters, 'IncidentEnergy'),     energy     = parameters.IncidentEnergy; end
  if isfield(parameters, 'IncidentWavevector'), wavevector = parameters.IncidentWavevector; end
  
  if isempty(lambda)
    if     ~isempty(energy)
      if ischar(energy), energy = get(s, energy); end
      lambda = sqrt(81.805./energy);
    elseif ~isempty(wavevector)
      if ischar(wavevector), wavevector = get(s, wavevector); end
      lambda = 2*pi./wavevector; 
    end
  end
  
  if isempty(lambda)
    disp([ mfilename ': WARNING: ' s.Tag ' ' s.Title ' undefined incident neutron wavelength/wavevector/energy.' ]);
    disp('    Define e.g. s.Wavelength=<lambda in Angs> or s.IncidentEnergy=<energy in meV>');
  end
  lambda = mean(lambda(:));
  if isempty(energy) || ischar(energy)     
    energy      = 81.805./lambda^2; 
    setalias(s, 'IncidentEnergy', energy, 'Incident Energy [meV]');
  end
  if isempty(wavevector)  || ischar(wavevector)   
    wavevector  = 2*pi./lambda;
    setalias(s, 'IncidentWavevector', energy, 'Incident Wavevector [Angs-1]');
  end
  
  % search for a sample-to-detector distance
  if isfield(parameters,'Distance'), 
    distance = mean(parameters.Distance(:));
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <sample-detector distance> =' num2str(distance) ' [m]' ]);
  end
  
  % search for the Channel Width
  if isfield(parameters, 'ChannelWidth'), 
    chwidth = mean(parameters.ChannelWidth(:));
    disp([ mfilename ': ' s.Tag ' ' s.Title ' using <channel width> =' num2str(chwidth) ' [time unit, s, ms or us]']);
  end
  

  
