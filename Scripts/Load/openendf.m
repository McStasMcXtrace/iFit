function out = openendf(filename)
%OPENENDF Open an Evaluated Nuclear Data File, display it
%        and set the 'ans' variable to an iData object with its content
% (c) E.Farhi, ILL. License: EUPL.

if ~isa(filename,'iData')
  out = iData(iLoad(filename,'ENDF'));  % no post-processing
else
  out = filename;
end
clear filename;

% must check if the input object should be split into temperatures
% the split of temperatures is done in openendf 
%  function t0=read_endf_mf7_array(t)
%  t0 = [];
%  disp(sprintf('%s: MF=%3i   MT=%i TSL T=%s', mfilename, t.MF, t.MT, mat2str(t.T)));
%  for index=1:numel(t.T)
%    t1     = t;
%    t1.T   = t.T(index);
%    t1.Title = [ t1.ZSYNAM ' T=' num2str(t1.T) ' [K] ' t1.description];
%    if t1.MT == 2 % Incoherent/Coherent Elastic Scattering
%      t1.NP  = t.NP(index);
%      t1.S   = t.S(index,:);
%      t1.INT = t.INT(index);
%      t1.Label = [ t1.ZSYNAM ' T=' num2str(t1.T) ' [K] TSL elastic' ];
%    elseif t1.MT == 4 % Incoherent Inelastic Scattering
%      t1.Sab = t.Sab(:,:,index);
%      t1.Teff     = t.Teff(index,:);
%      t1.Label = [ t1.ZSYNAM ' T=' num2str(t1.T) ' [K] TSL inelastic' ];
%    end
%    t0 = [ t0 t1 ];
%  end

if numel(out) > 1
  % handle input iData arrays
  for index=1:numel(out)
    out(index)  = feval(mfilename, out(index));
  end
  MF  = get(out,'MF');
  if numel(MF) == numel(out)
    out = out(find(MF ~= 1)); % not the MF1 section, which is stored into all objects
  end
elseif ~isempty(out) && isfield(out.Data,'MF') && isfield(out.Data,'MT')
  % set generic aliases
  setalias(out, 'MT',    'Data.MT',    'ENDF Section');
  setalias(out, 'MF',    'Data.MF',    'ENDF File');
  setalias(out, 'MAT',   'Data.MAT',   'ENDF Material number');
  setalias(out, 'Material','Data.ZSYNAM','ENDF Material description (ZSYNAM)');
  setalias(out, 'EDATE', 'Data.EDATE', 'ENDF Evaluation Date (EDATE)');
  setalias(out, 'charge','Data.ZA',    'ENDF material charge Z (ZA)');
  setalias(out, 'mass',  'Data.AWR',   'ENDF material mass A [g/mol] (AWR)');
  setalias(out, 'DescriptiveData', 'Data.DescriptiveData', 'ENDF DescriptiveData (MF1/MT451)');
  MT=out.Data.MT; MF=out.Data.MF;
  % assign axes: alpha, beta, Sab for MF7 MT4
  if MF == 1 && MT == 451
    setalias(out,'Header','Data.COMMENTS', 'Header from MF1/MT451');
  elseif MF == 7        % TSL:
    setalias(out,'weight','Data.AWR','Material mass [g/mol]');
    if MT == 2      % elastic
      if     out.Data.LTHR == 1
        setalias(out,'Signal','Data.S','S(E,T) Coherent Elastic Scattering Bragg Edges');
        setalias(out,'Energy','Data.E','Neutron Incident Energy [eV]');
        setalias(out,'Temperature','Data.T', 'Temperature [K]');
        setalias(out,'Scattering','coherent elastic');
        setaxis( out, 1, 'Energy');
      elseif out.Data.LTHR == 2
        setalias(out,'Signal','Data.W','Debye-Waller integral divided by the atomic mass [eV-1] ) (W)');
        setalias(out,'Sigma', 'Data.SB','Characteristic bound cross section [barns] (SB)');
        setalias(out,'Temperature','Data.T', 'Temperature [K]');
         setalias(out,'Scattering','incoherent elastic');
        setaxis( out, 1, 'Temperature');
      end
    elseif MT == 4  % inelastic
      setalias(out,'Signal','Data.Sab','S(alpha,beta,T) Inelastic Scattering');
      setalias(out,'alpha','Data.alpha','Unit-less wavevector [h2q2/AkT]');
      setalias(out,'beta', 'Data.beta', 'Unit-less energy [hw/kT]');
      setalias(out,'Temperature','Data.T', 'Temperature [K]');
      setalias(out,'classical','~this.Data.LASYM', 'classical/symmetric[1] or quantum/asymmetric[0]');
      setalias(out,'weight','this.Data.B(3)','Standard material mass [g/mol]');
      setalias(out,'sigma_free','this.Data.B(1)','total free atom cross section [barns]');
      setalias(out,'sigma_inc','this.Data.B(1)*(this.Data.B(3)+1)^2/this.Data.B(3)^2','total bound cross section [barns]');
      setalias(out,'multiplicity','this.Data.B(6)','the number of principal scattering atoms in the material');
      setalias(out,'B', 'Data.B', 'Analytical model physical constants (B)');
      setalias(out,'Scattering','(in)coherent inelastic');
      setaxis( out, 2, 'alpha');
      setaxis( out, 1, 'beta');
    end
  end
end

if ~nargout
  figure; subplot(log10(out));
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end
