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
  if ~isfield(out, 'MAT') && ~isempty(findfield(out, 'MAT','first case exact'))
    setalias(out, 'MAT',      findfield(out, 'MAT','first case'),   'ENDF Material number');
  end
  if ~isfield(out, 'ZSYMAM') && ~isempty(findfield(out, 'ZSYMAM','first case exact'))
    setalias(out, 'Material', findfield(out, 'ZSYMAM','first'),'ENDF Material description (ZSYMAM)');
  end
  if ~isfield(out, 'EDATE') && ~isempty(findfield(out, 'EDATE','first exact'))
    setalias(out, 'EDATE',    findfield(out, 'EDATE','first'), 'ENDF Evaluation Date (EDATE)');
  end
  if ~isfield(out, 'ZA') && ~isempty(findfield(out, 'ZA','first case exact'))
    setalias(out, 'charge',   findfield(out, 'ZA','first case'),    'ENDF material charge Z (ZA)');
  end
  if ~isfield(out, 'AWR') && ~isempty(findfield(out, 'AWR','first case exact'))
    setalias(out, 'mass',  findfield(out, 'AWR','first case'),   'ENDF material mass A [g/mol] (AWR)');
  end
  if ~isfield(out, 'DescriptiveData') && ~isempty(findfield(out, 'DescriptiveData','first case'))
    setalias(out, 'DescriptiveData', findfield(out, 'DescriptiveData','first case'), 'ENDF DescriptiveData (MF1/MT451)');
  end
  if isfield(out, 'info')
    setalias(out, 'DescriptiveData', findfield(out, 'info','first case'), 'ENDF DescriptiveData (MF1/MT451)');
  elseif ~isfield(out, 'description') && ~isempty(findfield(out, 'description','first case'))
    setalias(out, 'DescriptiveData', findfield(out, 'description','first case'), 'ENDF DescriptiveData (MF1/MT451)');
  end
  % search for aliases that we need
  tokens={'S','E','T','W','SB','Sab','alpha','beta','NP','INT','MT','MF'};
  for index=1:numel(tokens)
    tok = tokens{index};
    if ~isfield(out, tok) && ~isempty(findfield(out, tok,'first case exact'))
      setalias(out, tok, findfield(out, tok,'first case exact'));
    end
  end
  
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
  % must check if the input object should be split into temperatures
  if size(out, 3) > 1
    out = read_endf_mf7_array(out);
  end
elseif ~isempty(out) && isfield(out.Data, 'info')
  % we search for other structure members in Data and split the object into different sections
  out_array = [];
  f = fieldnames(out.Data);
  for index=1:numel(f)
    if strcmp(f{index}, 'info') continue; end
    section = out.Data.(f{index});
    if ~isstruct(section), continue; end
    % we detect some sections to keep
    if strcmp(f{index}, 'thermal_elastic')
      section.MF = 7; section.MT = 2;
    elseif strcmp(f{index},'thermal_inelastic')
      section.MF = 7; section.MT = 4;
    end
    info      = out.Data.info;
    out0      = copyobj(out); % out0.Data = []; 
    out0.Data = section;
    setalias(out0, 'info', info, 'ENDF General information (MF1)');
    out0      = feval(mfilename, out0);
    if numel(out0) > 1 || (~isempty(out0) && ~strcmp(getalias(out0,'Signal'), get(out0,'Signal')))
      out_array = [ out_array out0 ];
    end
    
  end
  out = out_array;
  
end % if not empty

if ~nargout
  figure; subplot(log10(out));
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

% ------------------------------------------------------------------------------
function t0=read_endf_mf7_array(t)
  % read_endf_mf7_array: split an array of ENDF entries vs Temperature
  t0 = [];
  disp(sprintf('%s: MF=%3i   MT=%i TSL T=%s', mfilename, t.MF, t.MT, mat2str(t.T)));
  for index=1:numel(t.T)
    t1     = copyobj(t);
    t1.T   = t.T(index);
    if isfield(t1, 'ZSYMAM')
      t1.Title = [ t1.Title ' ' t1.ZSYMAM ];
    end
    t1.Title = [ t1.Title ' T=' num2str(t1.T) ' [K]' ];
    if isfield(t1, 'description')
      t1.Title = [ t1.Title ' ' t1.description ];
    end
    if isfield(t1, 'ZSYMAM') t1.Label = t1.ZSYMAM; end
    t1.Label = [ t1.Label ' T=' num2str(t1.T) ' [K]' ];
    if t1.MT == 2 % Incoherent/Coherent Elastic Scattering
      t1.NP  = t.NP(index);
      t1.S   = t.S(index,:);
      t1.INT = t.INT(index);
      t1.Label = [ t1.Label ' TSL elastic' ];
      setaxis(t1, 3, t1.T);
    elseif t1.MT == 4 % Incoherent Inelastic Scattering
      if isfield(t1, 'Sab')
        t1.Sab = t.Sab(:,:,index);
      elseif isfield(t1, 'scattering_law')
        t1.Sab = t.scattering_law(:,:,index);
      end
      
      if isstruct(t.Teff) && isfield(t.Teff, 'y');
        t1.Teff     = t.Teff.y(index);
      else
        t1.Teff     = t.Teff(index);
      end
      t1.Label = [ t1.Label ' TSL inelastic' ];
      setaxis(t1, 3, t1.T);
    end
    t1 = iData(t1); % check axis and possibly transpose
    t0 = [ t0 t1 ];
  end
