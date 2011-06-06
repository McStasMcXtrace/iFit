function pars=iFuncs_private_guess(x, signal, parameter_names, dim)
% pars=iFuncs_private_guess(signal, parameter_names, dim)
%   guess private function to compute an estimate of parameters
%   given a signal and the model parameter names.
%   The following parameters can be set:
%     amplitude, intensity
%     peak width
%     peak center
%     signal background
%
% input:  a: signal array (double)
%         parameter_names: names of model parameters (cellstr)
%         dim: dimension to use. Default is 1 (int)
% output: pars: parameter values (double vector)
%
% Version: $Revision: 1.3 $

  pars=[];
  if nargin < 3, return; end
  if nargin < 4, dim=1;   end
  if isempty(parameter_names), return; end
  if isempty(x), x=1:length(signal); end
  if dim <= 0, dim=1; end

  pars            = zeros(size(parameter_names));
  parameter_names = lower(parameter_names);

  if ~isempty(signal)
    [sigma, position, amplitude, baseline] = iFuncs_private_findpeaks(signal, dim, 0);
  else
    % only the x axis is given: use default values from the axis extent
    sigma    = std(x)/2;
    position = round(length(x)/2);
    amplitude= 1;
    baseline = 0.1;
  end

  % sort peaks by amplitude
  [dummy,sorti] = sort(amplitude);
  sorti=sorti(end:-1:1);                  % descending amplitude
  amplitude = amplitude(sorti);
  sigma     = sigma(sorti);
  position  = position(sorti);

  % assign parameter guessed values according to their names
  for index=1:length(amplitude)
    % search for names that match a pattern, and not set previously
    for index_p=1:length(parameter_names)
      if pars(index_p) ~= 0, continue; end % was set before, go on with further parameters
      % test parameter names
      if     ~isempty(strfind(parameter_names{index_p}, 'amplitude')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'intensity')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'height')) ...
        pars(index_p) = amplitude(index);
      elseif ~isempty(strfind(parameter_names{index_p}, 'width')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'tau')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'damping'))
        pars(index_p) = sigma(index)*mean(diff(x))/2; 
      elseif ~isempty(strfind(parameter_names{index_p}, 'centre')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'center')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'position'))
        pars(index_p) = x(round(position(index))); 
      elseif ~isempty(strfind(parameter_names{index_p}, 'background')) ...
        |    ~isempty(strfind(parameter_names{index_p}, 'constant'))
        pars(index_p) = mean(baseline); 
      end
    end
  end
end % iFuncs_private_guess
