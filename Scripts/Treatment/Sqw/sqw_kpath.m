function S = sqw_kpath(f, qLim, E)
% sqw_kpath: evaluates a 4D S(q,w) model along specified k-path
%
%    sqw_kpath(f, kpath, w)
%    The k-path can be given as a cell containing 3-values (HKL) vectors, or
%      a n x 3 matrix, each row being a HKL location.
%    The energy range can be entered as a vector, or a [min max] pair.
%    Refer to https://en.wikipedia.org/wiki/Brillouin_zone and
%             http://lamp.tu-graz.ac.at/~hadley/ss1/bzones/ to get the standard 
%    points in the Brillouin zone.
%
% Example:
%   S = sqw_kpath(sqw_cubic_monoatomic, '', [0 10]); plot(log10(S));
%
% input:
%   f:    a 4D HKLE model S(q,w) (iFunc)
%   path: a list of k-locations given as a cell/array of HKL locations (cell or matrix)
%   w:    a vector of energies (4-th axis) for which to evaluate the model (double)
%
% output:
%   S:    the dispersion W(HKL) computed along the path (iData)
%
% Version: $Date$
% See also sqw_cubic_monoatomic, sqw_sine3d, sqw_vaks, sqw_spinw
%   sqw_powder, <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

  S = [];
  if nargin == 0, return; end
  
  if ndims(f) ~= 4 || ~isa(f,'iFunc')
    disp([ mfilename ': Invalid model dimension. Should be iFunc 4D. It is currently ' class(f) ' ' num2str(ndims(f)) 'D' ]);
    return
  end
  
  if nargin < 2, qLim = []; end
  if nargin < 3, E=[]; end
  if isempty(E)
    % make a quick evaluation in order to get the maxFreq
    qh=linspace(0.01,1.5,10);qk=qh; ql=qh; w=linspace(0.01,1000,100);
    [S,f] = feval(f, f.p, qh,qk,ql',w);
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
  if ischar(qLim), crystalsystem = qLim; qLim = {}; else crystalsystem = ''; end
  if isempty(crystalsystem) && isfield(f.UserData,'properties') ...
  && isfield(f.UserData.properties, 'spacegroup')
    spacegroup = regexp(f.UserData.properties.spacegroup,'\(([^:]*)\)','tokens');
    if ~isempty(spacegroup)
      spacegroup = str2double(spacegroup{1});
      if 195 <= spacegroup && spacegroup <= 230
        crystalsystem = 'Cubic';
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
      disp([ mfilename ': Using crystal system ' crystalsystem ', spacegroup ' f.UserData.properties.spacegroup ])
    end
  end
  
  % set default qLim path according to space group/crystal system
  if isempty(qLim)
    % we use when possible Path: Delta-Sigma-Lambda
    if isempty(crystalsystem), crystalsystem = 'orthorhombic'; end
    [qLim,xticks]=sqw_kpath_crystalsystem(crystalsystem);
  else xticks=[];
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
  % assemble all HKLw points
  h=[]; k=[]; l=[]; w=[]; index=1;
  for i=1:numel(H)
    for j=1:numel(E)
      h(index) = H(i);
      k(index) = K(i);
      l(index) = L(i);
      w(index) = E(j);
      index=index+1;
    end
  end
  % now we evaluate the model, without the intensities
  if isfield(f.UserData, 'properties') && isfield(f.UserData.properties, 'b_coh')
    b_coh = f.UserData.properties.b_coh;
  else b_coh = 0;
  end
  % unactivate intensity estimate, else some modes are not visible from |Q.e|
  f.UserData.properties.b_coh = 0; 
  [S, f, ax, name] = feval(f, f.p, h,k,l,E);
  f.UserData.properties.b_coh = b_coh;
  
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

  
  % create an iData
  x = linspace(0, numel(qLim)-1, numel(ax{1}));
  S = iData(E,x,S);
  % set title, labels, ...
  title(S, 'Model value along path');
  S.Title = [ f.Name ' along path' ];
  S=setalias(S, 'HKL',     qOut, 'HKL locations along path');
  S=setalias(S, 'kpoints', cell2mat(qLim'), 'set of HKL locations which form path');
  if ~isempty(xticks)
    S=setalias(S, 'XTickLabel', xticks, 'HKL labels');
  end
  
  xlabel(S,xlab);
  ylabel(S,'Energy');
  
  if ~isempty(inputname(1))
    assignin('caller',inputname(1),f); % update in original object
  end
  
% ----------------------------------------------------------------------------
function [qLim,lab]=sqw_kpath_crystalsystem(crystalsystem)
  % set a path depending on the crystal system
  qLim = {};
  switch lower(crystalsystem)
  case 'cubic'
    Gamma=[0,     0,     0    ];
    X=    [0,     0 / 2, 1 / 2];
    R=    [1 / 2, 1 / 2, 1 / 2];
    M=    [0 / 2, 1 / 2, 1 / 2];
    qLim = {Gamma X M Gamma R M};
    lab  ='Gamma X M Gamma R M';
  case 'fcc'
    Gamma=[0,     0,     0    ];
    X=    [1 / 2, 0,     1 / 2];
    W=    [1 / 2, 1 / 4, 3 / 4];
    K=    [3 / 8, 3 / 8, 3 / 4];
    U=    [5 / 8, 1 / 4, 5 / 8];
    L=    [1 / 2, 1 / 2, 1 / 2];
    qLim = { Gamma X W K Gamma L U W };
    lab  = 'Gamma X W K Gamma L U W';
  case 'bcc'
    Gamma=[0,      0,     0    ];
    H=    [1 / 2, -1 / 2, 1 / 2];
    N=    [0,      0,     1 / 2];
    P=    [1 / 4,  1 / 4, 1 / 4];
    qLim = { Gamma H N Gamma P H };
    lab  = 'Gamma H N Gamma P H';
  case 'hexagonal'
    Gamma= [0,      0,       0   ];
    M=     [0,      1 / 2,   0   ];
    K=     [-1 / 3, 1 / 3,   0   ];
    A=     [0,      0,     1 / 2 ];
    L=     [0,     1 / 2,  1 / 2 ];
    H=     [-1 / 3, 1 / 3, 1 / 2 ];
    qLim = { Gamma A L M Gamma K H };
    lab  = 'Gamma A L M Gamma K H';
  case 'tetragonal'
    Gamma= [0,      0,       0   ];
    X=     [1 / 2,  0,       0   ];
    M=     [1 / 2,  1 / 2,   0   ];
    Z=     [0,      0,     1 / 2 ];
    R=     [1 / 2,  0,     1 / 2 ];
    A=     [1 / 2,  1 / 2, 1 / 2 ];
    qLim = { Gamma X M Gamma Z R A };
    lab  = 'Gamma X M Gamma Z R A';
  case 'orthorhombic'
    Gamma= [0,      0,       0   ];
    R=     [1 / 2,  1 / 2, 1 / 2 ];
    S=     [1 / 2,  1 / 2,   0   ];
    T=     [0,      1 / 2, 1 / 2 ];
    U=     [1 / 2,  0,     1 / 2 ];
    X=     [1 / 2,  0,       0   ];
    Y=     [0,      1 / 2,   0   ];
    Z=     [0,      0,     1 / 2 ];
    qLim = { Gamma Y S X Gamma Z T R U };
    lab  = 'Gamma Y S X Gamma Z T R U';
  end

