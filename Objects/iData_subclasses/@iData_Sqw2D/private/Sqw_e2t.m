function s=Sqw_e2t(s)
% convert S(xx,w) to S(xx,t). Requires wavelength, chwidth, distance

  if isempty(s), return; end
  [s,lambda,distance,chwidth] = Sqw_search_lambda(s);

  disp([ mfilename ': ' s.Tag ' ' s.Title ' Converting Axis 1 "' ...
    label(s,1) '": energy [meV] to time [sec].' ]);
 
  if isempty(lambda),   lambda   = 2.36; end
  if isempty(distance), distance = 4; end
  if isempty(chwidth),  chwidth  = 4; end
  

  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
  Ki   = 2*pi/lambda;
  Vi   = K2V*Ki;
  Ei   = VS2E*Vi.^2;
  
  % compute final energy
  hw   = getaxis(s, 1);
  Ef   = Ei - hw;
  t_sample_detector = distance/Vi;  % this is the time for the elastic peak
  
  t  = t_sample_detector./sqrt(Ef./Ei);
  t  = t / chwidth;
  
  dtdE    = t./(2.*Ef)*1e6;   % to be consistent with vnorm abs. calc.
  kikf    = sqrt(Ei./Ef);     % all times above calculated in sec.
  dtdEkikf= dtdE.*kikf;
  s       = s./dtdEkikf/chwidth;
  
  setaxis(s, 1, t);
  label(s, 1, 'Time channel [channel]');
