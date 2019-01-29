function s4 = iFunc_Sqw2Dto4D(s, B)
% iFunc_Sqw2Dto4D: convert a Sqw2D model to a 4D one, given the lattice parameters

% B not given: search for lattice parameters/B matrix in object
if nargin < 2 
  UD = s.UserData;
  if isfield(UD, 'reciprocal_cell')
    B = UD.reciprocal_cell;
  elseif isfield(UD, 'properties') && isfield(UD.properties, 'reciprocal_cell')
    B = UD.properties.reciprocal_cell;
  elseif ~isempty(findfield(s, 'reciprocal_cell'))
    index = findfield(s, 'reciprocal_cell','first cache exact');
    if iscell(index), index=index{1}; end
    B = get(sa, index);
  elseif isfield(UD, 'cell')
    B = UD.reciprocal_cell;
  elseif isfield(UD, 'properties') && isfield(UD.properties, 'cell')
    B = UD.properties.reciprocal_cell;
  elseif ~isempty(findfield(s, 'cell'))
    index = findfield(s, 'cell','first cache exact');
    if iscell(index), index=index{1}; end
    B = get(sa, index);
  else
    warning([ mfilename ': WARNING: no reciprocal_cell information found. Assuming cubic a=2*pi.' ]);
    B = eye(3); % assume cubic, a=b=c=2*pi, 90 deg, then a*=2pi/a=1...
  end
end

% compute B matrix if needed
if isvector(B) && numel(B) == 6
  % compute B matrix from [a b c alpha beta gamma]
  s.UserData.cell = B;
  alpha=s.UserData.cell(4);
  beta =s.UserData.cell(5);
  gamma=s.UserData.cell(6);
  a_vec=s.UserData.cell(1)*[1; 0; 0];
  b_vec=s.UserData.cell(2)*[cosd(gamma); sind(gamma); 0];
  c1=cosd(beta); 
  c2=(cosd(alpha)-cosd(gamma)*cosd(beta))/sind(gamma); 
  c3=sqrt(1-c1^2-c2^2);
  c_vec=s.UserData.cell(3)*[c1; c2; c3;];
  V=dot(a_vec,cross(b_vec,c_vec));
  B=2*pi*[cross(b_vec,c_vec) cross(c_vec,a_vec) cross(a_vec,b_vec)]/V; % reciprocal basis, as columns
  s.UserData.volume = V;
elseif numel(B) ~= 9 || any(size(B) ~= 3)
  warning([ mfilename ': Conversion from 2D model to 4D[rlu] requires reciprocal/lattice cell parameters. Use: iData_Sqw4D(iData_Sqw2D, B or [a b c alpha beta gamma]). Assuming cubic.' ]);
  warning([ mfilename ': WARNING: no reciprocal_cell information found. Assuming cubic a=2*pi.' ]);
  B = eye(3); % assume cubic, a=b=c=2*pi, 90 deg, then a*=2pi/a=1...
  s.UserData.cell = [ 2*pi 2*pi 2*pi 90 90 90  ];
end
s.UserData.reciprocal_cell = B;
s.UserData.B = B;

% now create new object from 2D one
% prepended code: get [x,y,z] = [hkl], compute [qx,qy,qz], then |q| in Angs-1

% the template below creates HKL = [ x(:) y(:) z(:) ] in [rlu]
script_hkl = sqw_phonons_templates;

s4 = copyobj(iFunc(s)); % copy all internal stuff, including UserData
s4.Dimension  = 4;
e=s.Expression;
s4.Expression = { ...
  script_hkl{:} ...
  'B = this.UserData.reciprocal_cell;', ...
  'qxyz = B*HKL'';' ...
  'qx=qxyz(1,:); qy=qxyz(2,:); qz=qxyz(3,:);' ...
  'q=sqrt(qx.^2+qy.^2+qz.^2); % norm(q)' ...
  'x0=x; y0=y; signal=0;' ...
  'x =q(:); y =t(:);' ...
  'if isvector(x) && isvector(y) && numel(x) ~= numel(y), [y,x] = meshgrid(y,x); end' ...
  '% now evaluate the 2D model S(q,w)' ...
  e{:} ...
  'signal(~isfinite(signal) | signal < 0 | signal > 1e10) = 0;' ...
  'if ~isempty(resize_me) && prod(resize_me) == numel(signal)' ...
  'signal = reshape(signal, resize_me); % initial 4D cube dimension = [ size(x) numel(t) ]' ...
  'end' ...
  'x=x0; y=y0;' ...
};


