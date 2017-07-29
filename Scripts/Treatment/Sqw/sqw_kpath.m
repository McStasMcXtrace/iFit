function [S, qLim, fig] = sqw_kpath(f, qLim, E, options)
% sqw_kpath: evaluates a 4D S(q,w) model along specified k-path / bands
%
%    sqw_kpath(f, kpath, w, options)
%      The k-path can be given as a cell containing 3-values (HKL) vectors, or
%        a n x 3 matrix, each row being a HKL location.
%      The energy range can be entered as a vector, or a [min max] pair.
%      Refer to https://en.wikipedia.org/wiki/Brillouin_zone and
%             http://lamp.tu-graz.ac.at/~hadley/ss1/bzones/ to get the standard 
%      points in the Brillouin zone.
%
%      The bands intensity is computed along the path. For neutron scattering,
%        as all standard points in the Brillouin zone are centered around 0 and
%        the intensity is proportional to |Q.e|^2, all transverse modes are extinct.
%      To get the proper band intensity around a given Bragg peak G, the k-path
%        should be shifted by G.
%
%    sqw_kpath(f);
%      plots the dispersion curves following the BZ points for the 
%      crystal spacegroup, and using the maximum excitation energy.
%
%    [S,k] = sqw_kpath(...)
%      returns the dispersion curves data set (iData), and the k-path used.
%
%    When the command is followed by ';' or options contains 'plot', a plot is 
%    generated.
%
% Example:
%   S = sqw_kpath(sqw_cubic_monoatomic, '', [0 10]); plot(log10(S));
%   S = sqw_phonons('POSCAR_Al','emt','metal'); sqw_kpath(S);
%
% input:
%   f:    a 4D HKLE model S(q,w) (iFunc)
%   path: a list of k-locations given as a cell/array of HKL locations (cell or matrix)
%         can also be 'Cubic','Hexagonal','Trigonal','Tetragonal','Orthorhombic',
%                     'Monoclinic','Triclinic','fcc','bcc'
%         can also be given as a list of BZ points, such as:
%           {'Gamma' 'K' 'M' 'Gamma' 'A' 'H' 'L' 'A' }
%           which are then defined according to the space group, e.g:
%             f.UserData.properties.spacegroup_number
%         can be given as a list of HKL locations, such as {[0,0,0],[1,1,1],[1,1,0]}
%           or as an array e.g. [0 0 0 ; 1 1 1 ; 1 1 0 ; 0 0 0 ; 0 0 1]
%         The path can also be given as a Bragg peak location, e.g. [0 0 1] to compute
%           the dispersion around that location using the default K-path in the BZ.
%         when not given or empty, this is guessed from the crystal spacegroup.
%
%   w:    a vector of energies (4-th axis) for which to evaluate the model (double)
%         can also be given as Emax, or [ Emin Emax ]
%         when not given or empty, the maximum excitation energy is used.
%   options: plotting option (string). Can contain 'plot' 'THz','cm-1','meV'.
%         The 'plot' option also displays the density of states, when available.
%
% output:
%   S:    the dispersion W(HKL) computed along the path (iData)
%   k:    the k-path used to generate the dispersion curves (matrix)
%   fig:  figure handle generated when options contains 'plot' or no output.
%
% Version: $Date$
% See also sqw_cubic_monoatomic, sqw_sine3d, sqw_vaks, sqw_spinw
%   sqw_powder, <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

  S = []; fig = [];
  if nargin == 0, return; end
  if nargin < 2, qLim = []; end
  if nargin < 3, E=[]; end
  if nargin < 4, options=[]; end
  
  % handle input array
  if numel(f) > 1
    for index=1:numel(f)
      S = [ S sqw_kpath(f(index), qLim, E, options) ];
    end
    return
  end
  
  if ndims(f) ~= 4 || ~isa(f,'iFunc')
    disp([ mfilename ': Invalid model dimension. Should be iFunc 4D. It is currently ' class(f) ' ' num2str(ndims(f)) 'D' ]);
    return
  end
  
  % default options in qLim or E ?
  if ischar(qLim) && any(~cellfun(@isempty, ...
     strfind(lower({'plot' 'THz','cm-1','meV','cm','hz'}),strtok(lower(qLim)))))
    options = qLim; qLim = '';
  elseif ischar(E) && any(~cellfun(@isempty, ...
     strfind(lower({'plot' 'THz','cm-1','meV','cm','hz'}), strtok(lower(E)))))
    options = E; E = [];
  end
  
  if isempty(E)
    % make a quick evaluation in order to get the maxFreq
    qh=linspace(0.01,1.5,10);qk=qh; ql=qh; w=linspace(0.01,1000,100);
    % this also computes the DOS when not there yet. In case it was not there before, we remove it.
    is_dos_there = isfield(f.UserData,'DOS') && ~isempty(f.UserData.DOS) ...
      && isa(f.UserData.DOS, 'iData') && prod(size(f.UserData.DOS)) > 10000;
    [S,f] = feval(f, f.p, qh,qk,ql',w);
    if ~is_dos_there, f.UserData.DOS=[]; end
    % search for maxFreq if exists
    if isfield(f.UserData, 'maxFreq')
      E = max(f.UserData.maxFreq)*1.2;
    else 
      S = iData(w,squeeze(sum(sum(sum(S)))));
      [half_width, center] = std(S);
      E=center+half_width;
    end
  end
  if numel(E) == 1
    E = linspace(0.01,E, 100);
  elseif numel(E) == 2
    E = linspace(min(E),max(E), 100);
  end

  % look for space group when available
  spacegroup = '';
  if ischar(qLim), crystalsystem = qLim; qLim = {}; else crystalsystem = ''; end
  if isempty(crystalsystem) && isfield(f.UserData,'properties') ...
  && isfield(f.UserData.properties, 'spacegroup')
    spacegroup = regexp(f.UserData.properties.spacegroup,'\(([^:]*)\)','tokens');
    if iscell(spacegroup) && ~isempty(spacegroup), spacegroup = spacegroup{1}; end
    if isempty(spacegroup) && isfield(f.UserData.properties, 'spacegroup_number')
      spacegroup = f.UserData.properties.spacegroup_number;
    end
    if ~isempty(spacegroup)
      if iscell(spacegroup), spacegroup = spacegroup{1}; end
      if ischar(spacegroup), spacegroup = str2double(spacegroup); end
      if 195 <= spacegroup && spacegroup <= 230
        if lower(f.UserData.properties.spacegroup(1)) == 'f'
          crystalsystem = 'fcc';
        elseif lower(f.UserData.properties.spacegroup(1)) == 'i'
          crystalsystem = 'bcc';
        else
          crystalsystem = 'Cubic';
        end
      elseif 168 <= spacegroup && spacegroup <= 194
        crystalsystem = 'Hexagonal';
      elseif 143 <= spacegroup && spacegroup <= 167
        crystalsystem = 'Trigonal';
      elseif 75 <= spacegroup && spacegroup <= 142
        crystalsystem = 'Tetragonal';
      elseif 16 <= spacegroup && spacegroup <= 74
        crystalsystem = 'Orthorhombic';
      elseif 3 <= spacegroup && spacegroup <= 15
        crystalsystem = 'Monoclinic';
      elseif 1 <= spacegroup && spacegroup <= 2
        crystalsystem = 'Triclinic';
      end
    end
    if ~isempty(crystalsystem)
      disp([ mfilename ': Using crystal system ' crystalsystem ', spacegroup ' f.UserData.properties.spacegroup ' (' num2str(spacegroup) ')' ])
    end
  end
  
  % define default path and BZ points
  if isempty(crystalsystem), crystalsystem = 'orthorhombic'; end
  [qLim0,lab, points]=sqw_kpath_crystalsystem(crystalsystem);
  
  % set default qLim path according to space group/crystal system
  xticks=[];

  if (numel(qLim) == 1 || numel(qLim) == 3) && isnumeric(qLim)
    if numel(qLim) == 1
      qLim = [ qLim qLim qLim ];
    end
    for i=1:numel(qLim0); 
      qLim0{i} = qLim0{i}+qLim; 
    end
    qLim = qLim0;
    xticks = lab;
  elseif isempty(qLim)
    % we use when possible Path: Delta-Sigma-Lambda
    qLim = qLim0;
    xticks = lab;
  elseif iscellstr(qLim) % a list of BZ point names
    for index=1:numel(qLim)
      if isfield(points, qLim{index})
        
        xticks = [ xticks qLim{index} ' ' ];
        qLim{index} = points.(qLim{index});
      else
        disp([ mfilename ': Invalid BZ point ' qLim{index} ' in specified path. Removing it.' ]);
        qLim{index}=[];
      end
    end
    qLim = qLim(~cellfun(@isempty, qLim));
  end
  if ~iscell(qLim) && ndims(qLim) == 2 && isnumeric(qLim)
    if size(qLim,1) == 3 && size(qLim,2) ~= 3
      qLim = Lim';
    end
    if size(qLim,2) == 3
      q = cell(1,size(qLim,1));
      for index=1:size(qLim,1)
        q{index} = qLim(index,:);
      end
      qLim = q;
    end
  end
  if ~iscell(qLim) || isempty(qLim)
    disp([ mfilename ': Invalid kpath argument. Should be a cell or nx3 matrix with HKL locations. It is currently ' class(qLim) ]);
    return
  end
  
  % we use SpinW qscan function to generate the k-points along given directions
  qOut = sw_qscan(qLim);
  
  % now we compute the model value along the path
  H = qOut(1,:);
  K = qOut(2,:);
  L = qOut(3,:);

  % assemble all HKLw points (event list)
  h=ones(numel(E),1)*H; h=h(:)';
  k=ones(numel(E),1)*K; k=k(:)';
  l=ones(numel(E),1)*L; l=l(:)';
  w=E'*ones(1,numel(H)); w=w(:)';

  % now we evaluate the model, without the intensities
  if isfield(f.UserData, 'properties') && isfield(f.UserData.properties, 'b_coh')
    b_coh = f.UserData.properties.b_coh;
  else b_coh = 0;
  end
  % unactivate intensity estimate, else some modes may not be visible from |Q.e|
  if isfield(f.UserData,'properties')
    f.UserData.properties.b_coh = 0; 
  end

  [S, f, ax, name] = feval(f, f.p, h,k,l,E);
  if isfield(f.UserData,'properties')
    f.UserData.properties.b_coh = b_coh;
  end
  
  if isfield(f.UserData,'properties') && isfield(f.UserData.properties,'chemical_formula')
    chem = [ ' ' f.UserData.properties.chemical_formula ];
  else chem = ''; end
  
  % retain only HKL locations (get rid of optional n)
  if isscalar(qLim{end}), 
    n = qLim{end}; 
    qLim(end) = [];
  end
  
  xlab = '';
  for index=1:numel(qLim)
    if index==1, xlab = sprintf('[%g %g %g]', qLim{index});
    else         xlab = [xlab sprintf(' : [%g %g %g]', qLim{index})];
    end
  end
  
  % now we generate a 2D iData
  % S = reshape(S, [  numel(ax{1}) numel(E) ]);
  
  if strfind(lower(options), 'hz')
    factor = 0.2418;  % meV -> Thz
    unit = 'THz';
  elseif strfind(lower(options), 'cm')
    factor = 8.0657;  % meV -> cm-1
    unit = 'cm-1';
  else
    factor = 1;
    unit = 'meV';
  end
  
  % create an iData
  x = linspace(0, numel(qLim)-1, numel(ax{1}));
  S = iData(E*factor,x,S);
  % set title, labels, ...
  if ~isempty(spacegroup)
    spacegroup = sprintf(' (sg %i)', spacegroup);
  end
  title(S, [ 'Model value along path' chem spacegroup ]);
  S.Title = [ f.Name ' along path. ' crystalsystem chem spacegroup ];
  S=setalias(S, 'HKL',     qOut, 'HKL locations along path');
  S=setalias(S, 'kpoints', cell2mat(qLim'), 'set of HKL locations which form path');
  if ~isempty(xticks)
    S=setalias(S, 'XTickLabel', xticks, 'HKL labels');
  end
  xlabel(S,xlab);
  ylabel(S,[ 'Energy ' unit ]);
  
  qLim = cell2mat(qLim');
  
  % plot results when no output

  if nargout == 0 || ~isempty(strfind(options, 'plot'))
    fig = figure; 
    if isfinite(max(S)) && max(S), plot(log10(S/max(S)),'view2'); 
    else plot(log10(S),'view2'); end
    axis tight
    add_contextmenu(gca)
    hold on
    if isfield(f.UserData,'FREQ')
      FREQ = f.UserData.FREQ*factor;
      if ~isempty(FREQ), plot(x, FREQ); end
    end
    % evaluate DOS if not done yet
    if ~isfield(f.UserData,'DOS') || isempty(f.UserData.DOS)
      qh=linspace(-.5,.5,30);qk=qh; ql=qh; w=linspace(0.01,50,11);
      F=iData(f,[],qh,qk,ql',w);
      clear F
    end
    if isfield(f.UserData,'DOS') && ~isempty(f.UserData.DOS)
      figure(fig);
      DOS = f.UserData.DOS;
      DOS{1} = DOS{1}*factor; % change energy unit
      xlabel(DOS,[ 'Energy ' unit ]);
      % rescale dispersion curves and use same limits for DOS
      p = get(gca,'Position'); p(3) = 0.6; set(gca,'Position', p);
      y = ylim(gca);
      a = axes('position', [ 0.8 p(2) 0.15 p(4) ]);
      % plot any partials first
      if isfield(f.UserData,'DOS_partials') && numel(f.UserData.DOS_partials) > 0
        d=f.UserData.DOS_partials;
        for index=1:numel(d)
          this_pDOS=d(index);
          this_pDOS{1} = this_pDOS{1}*factor;
          d(index) = this_pDOS;
        end
        h=plot(d);
        if iscell(h)
            if isnumeric(h{1}), h=cell2mat(h); 
            else h = [ h{:} ]; end
        end
        set(h,'LineStyle','--');
        hold on
      end
      % plot total DOS and rotate
      h=plot(DOS); set(h,'LineWidth',2);
      xlabel(''); ylabel('DOS'); title('');
      xlim(y);
      view([90 -90]);
    end
    hold off
  else fig = [];
  end
  
  if ~isempty(inputname(1))
    assignin('caller',inputname(1),f); % update in original object
  end
  
% ----------------------------------------------------------------------------
function [qLim,lab, points]=sqw_kpath_crystalsystem(crystalsystem)
  % set a path depending on the crystal system
  % may be given from spacegroup, or as a string:
  %   'Cubic','Hexagonal','Trigonal','Tetragonal','Orthorhombic','Monoclinic';
  %   'Triclinic','fcc','bcc'
  qLim = {};
  switch lower(crystalsystem)
  case 'cubic'
    Gamma=[0,     0,     0    ];
    X=    [0,     0 / 2, 1 / 2];
    R=    [1 / 2, 1 / 2, 1 / 2];
    M=    [0 / 2, 1 / 2, 1 / 2];
    qLim = {Gamma X M Gamma R M};
    lab  = 'Gamma X M Gamma R M';
  case 'fcc'
    Gamma=[0,     0,     0    ];
    X=    [1 / 2, 0,     1 / 2];
    W=    [1 / 2, 1 / 4, 3 / 4];
    K=    [3 / 8, 3 / 8, 3 / 4];
    U=    [5 / 8, 1 / 4, 5 / 8];
    L=    [1 / 2, 1 / 2, 1 / 2];
    qLim = { Gamma X U K Gamma L W X };
    lab  =  'Gamma X U K Gamma L W X';
  case 'bcc'
    Gamma=[0,      0,     0    ];
    H=    [1 / 2, -1 / 2, 1 / 2];
    N=    [0,      0,     1 / 2];
    P=    [1 / 4,  1 / 4, 1 / 4];
    qLim = { Gamma H N Gamma P H };
    lab  =  'Gamma H N Gamma P H';
  case {'hexagonal','trigonal'}
    Gamma= [0,      0,       0   ];
    M=     [0,      1 / 2,   0   ];
    K=     [-1 / 3, 1 / 3,   0   ];
    A=     [0,      0,     1 / 2 ];
    L=     [0,     1 / 2,  1 / 2 ];
    H=     [-1 / 3, 1 / 3, 1 / 2 ];
    qLim = { Gamma A L M Gamma K H };
    lab  =  'Gamma A L M Gamma K H';
  case 'tetragonal'
    Gamma= [0,      0,       0   ];
    X=     [1 / 2,  0,       0   ];
    M=     [1 / 2,  1 / 2,   0   ];
    Z=     [0,      0,     1 / 2 ];
    R=     [1 / 2,  0,     1 / 2 ];
    A=     [1 / 2,  1 / 2, 1 / 2 ];
    qLim = { Gamma X M Gamma Z R A };
    lab  =  'Gamma X M Gamma Z R A';
  otherwise % 'orthorhombic','triclinic','monoclinic'
    Gamma= [0,      0,       0   ];
    R=     [1 / 2,  1 / 2, 1 / 2 ];
    S=     [1 / 2,  1 / 2,   0   ];
    T=     [0,      1 / 2, 1 / 2 ];
    U=     [1 / 2,  0,     1 / 2 ];
    X=     [1 / 2,  0,       0   ];
    Y=     [0,      1 / 2,   0   ];
    Z=     [0,      0,     1 / 2 ];
    qLim = { Gamma Y S X Gamma Z T R U };
    lab  =  'Gamma Y S X Gamma Z T R U';
  end
  
  l    = textscan(lab,'%s','Delimiter', ' '); l = l{1};
  for f = l'
    points.(f{1}) = eval(f{1});
  end

function add_contextmenu(a)
  
  uicm = uicontextmenu;
  % menu Duplicate (axis frame/window)
  uimenu(uicm, 'Label', 'Duplicate View...', 'Callback', ...
     [ 'tmp_cb.g=gca;' ...
       'tmp_cb.f=figure; tmp_cb.c=copyobj(tmp_cb.g,gcf); ' ...
       'set(tmp_cb.c,''position'',[ 0.1 0.1 0.85 0.8]);' ...
       'set(gca,''XTickLabelMode'',''auto'',''XTickMode'',''auto'');' ...
       'set(gca,''YTickLabelMode'',''auto'',''YTickMode'',''auto'');' ...
       'set(gca,''ZTickLabelMode'',''auto'',''ZTickMode'',''auto'');']);
       
  uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');

  uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback', [ ...
    '[tmp_a,tmp_e]=view; if (tmp_a==0 & tmp_e==90) view(3); else view(2); end;' ...
    'clear tmp_a tmp_e; lighting none;alpha(1);shading flat;rotate3d off;axis tight;' ]);
  uimenu(uicm, 'Label','Smooth View','Callback', 'shading interp;');
  uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
  uimenu(uicm, 'Label','Transparency','Callback', 'alpha(0.7);');
  uimenu(uicm, 'Label', 'Linear/Log Intensity','Callback', 'if strcmp(get(gca,''zscale''),''linear'')  set(gca,''zscale'',''log''); else set(gca,''zscale'',''linear''); end');
  uimenu(uicm, 'Label','Linear/Log X axis', ...
  'Callback', 'if strcmp(get(gca,''xscale''),''linear'')  set(gca,''xscale'',''log''); else set(gca,''xscale'',''linear''); end');
  uimenu(uicm, 'Label','Linear/Log Y axis', ...
  'Callback', 'if strcmp(get(gca,''yscale''),''linear'')  set(gca,''yscale'',''log''); else set(gca,''yscale'',''linear''); end');
  uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
  

  uimenu(uicm, 'Separator','on','Label', 'About iFit/iData', ...
    'Callback',[ 'msgbox(''' version(iData,2) sprintf('. Visit <http://ifit.mccode.org>') ''',''About iFit'',''help'')' ]);
  set(a, 'UIContextMenu', uicm);
