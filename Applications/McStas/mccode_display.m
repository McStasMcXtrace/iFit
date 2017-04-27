function comps=mccode_display(model)
% runs the model in --trace mode and capture the output
% grab all MCDISPLAY lines

  % switch to trace mode
  model.UserData.options.trace  = 1;
  model.UserData.options.ncount = 1e3;
  % execute and get output
  output = evalc('val=model([],nan);');
  % check that we have MCDISPLAY words in there.
  if isempty(strfind(output, 'MCDISPLAY: '))
    disp([ mfilename ': Can not get MCDISPLAY tokens from the simulation. Something probably went wrong.' ]);
    return
  end

  output = textscan(output, '%s','Delimiter','\n\r');
  output = output{1};

  % first extract the portion 'start':'end'
  index_start = find(~cellfun(@isempty, strfind(output, 'MCDISPLAY: start')));
  index_end   = find(~cellfun(@isempty, strfind(output, 'MCDISPLAY: end')));
  if numel(index_start) ~= 1 || numel(index_end) ~= 1
    disp([ mfilename ': The MCDISPLAY section is invalid (incomplete or multiple). Aborting.' ]);
    return
  end

  % initiate component structures
  % we build a struct array, with one entry per component, and fields:
  %   name: name
  %   pos: pos(3)
  %   rot: rot(3,3)
  %   x,y,z: set of points
  comps = mcdisplay_get_components(output(1:index_start));
  % restrict output to the MCDISPLAY so that all searches are faster
  output = output(index_start:index_end); 

  % get the components in order, and identify the output section/lines
  % which are separated by e.g. MCDISPLAY: component <blah>
  index_mcdisplay_comp = find(~cellfun(@isempty, strfind(output, 'MCDISPLAY: component ')));
  if numel(index_mcdisplay_comp) ~= numel(comps)
    disp([ mfilename ...
      ': WARNING: not the same number of declared components (' num2str(numel(comps)) ...
      ') and MCDISPLAY sections ' num2str(numel(index_mcdisplay_comp)) ])
  end
  
  % extract the multiline and circle stuff in each component mcdisplay section
  for index=1:numel(index_mcdisplay_comp)
    if index < numel(index_mcdisplay_comp), 
      next = index_mcdisplay_comp(index+1);
    else 
      next = numel(output); end
    % get the MCDISPLAY section for a single component
    section = output(index_mcdisplay_comp(index):next);
    % then we get the multiline and circle commands in this section
    for token = {'multiline' ,'circle'}
      [x,y,z] = mcdisplay_get_token(section, token{1});
      comps(index).x = [ comps(index).x nan x ];
      comps(index).y = [ comps(index).y nan y ];
      comps(index).z = [ comps(index).z nan z ];
    end
  end
  
  clear output
  
  % transform the points and plot them
  f = figure('Name',[ 'Instrument: ' model.name ]);
  colors='bgrcmk';
  for index=1:numel(comps)
    comp = comps(index);
    disp([' Component: ' comp.name ' [' num2str(index) ']' ])
    r = [ comp.x ; comp.y ; comp.z ];
    if all(isnan(r)), continue; end
    R = comp.rot*r;
    x = R(1,:)+comp.pos(1);
    y = R(2,:)+comp.pos(2);
    z = R(3,:)+comp.pos(3);
    c = mod(comp.index, numel(colors)); c=colors(c+1);
    h = plot3(z,x,y, [ c '-' ]);
    popup=uicontextmenu;
    uimenu(popup,'label', comp.name,'ForeGroundColor',c);
    uimenu(popup,'label', [ 'AT: ' mat2str(comp.pos) ]);
    set(h,'uicontextmenu',popup);
    hold on
  end
  xlabel('Z [m]');
  ylabel('X [m]');
  zlabel('Y [m]')
  
  mp    = model.ParameterValues;
  names = model.Parameters;
  t = [];
  for index=1:numel(mp)
    if numel(mp) < index, val = []; else val = mp(index); end
    t = [ t names{index} '=' num2str(val) ];
  end
  title({ [ 'Instrument: ' model.name ] ; t });

end % mcdis





% ------------------------------------------------------------------------------
% initialize the component structures by searching name and pos/rot, e.g.
%   COMPONENT: "collimador_radial"
%   POS: 1.82089, 0, 19.6314, 0.822317, 0, -0.569029, -0, 1, 0, 0.569029, -0, 0.822317
function comps = mcdisplay_get_components(output)
  token = 'COMPONENT: "';
  index_token = find(~cellfun(@isempty, strfind(output, token)));
  comps = [];
  for index = index_token'
    this_line = output{index};
    compname  = strtok(this_line(numel(token):end),'"');
    pos = [];
    if index < numel(output)
      next_line = output{index+1};
      if strncmp(next_line, 'POS: ', 5), pos = str2num(next_line(6:end)); end
    end

    if ~isempty(compname) && numel(pos) == 12
      comp.name = compname;
      comp.pos  = pos(1:3);                 % absolute position
      comp.rot  = reshape(pos(4:end),3,3);  % absolute rotation matrix
      comp.index= numel(comps)+1;
      comp.x=[]; comp.y=[]; comp.z=[];
      comps = [ comps comp ];
    end
  end
end

% search for a plot command: multiline or circle. Get the local coords back.
function [X,Y,Z] = mcdisplay_get_token(output, token)
  index_token = find(~cellfun(@isempty, strfind(output, [ 'MCDISPLAY: ' token ])));
  X=[]; Y=[]; Z=[];
  if isempty(index_token), return; end
  for index = index_token
    % in each multiline, we replace the search token and execute the remaining part in eval
    this_line = output{index};
    this_line = strrep(this_line, 'MCDISPLAY: ','[x,y,z]=');
    eval([ this_line ';' ]);  % executes multiline(..), return a set of points (local coords)
    X = [ X x ];
    Y = [ Y y ];
    Z = [ Z z ];
  end
end

% get the multiline and circle arguments back
function [X,Y,Z]=multiline(npoints, varargin)
  X=[]; Y=[]; Z=[];
  for index=1:3:numel(varargin)
    X=[ X varargin{index} ];
    Y=[ Y varargin{index+1} ];
    Z=[ Z varargin{index+2} ];
  end
end

function [X,Y,Z]=circle(plane, x0,y0,z0, radius)
  % we create a set of points
  if radius ~=0
    phi=linspace(0,2*pi, 36); % 36 points along the circle
  else phi=0; end
  x=radius*sin(phi);
  y=radius*cos(phi);
  zero = 0*x;
  switch plane
  case {'xy','yx'}
    X=x; Y=y; Z=zero;
  case {'xz','zx'}
    X=x; Y=zero; Z=y;
  case {'zy','yz'}
    X=zero; Y=x; Z=y;
  otherwise
    X=[]; Y=[]; Z=[]; 
    disp([ mfilename ': unknown plane: ' plane ]);
    return
  end
  X=X+x0; Y=Y+y0; Z=Z+z0;
end

