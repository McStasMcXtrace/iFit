function filename = iData_private_saveas_x3d(a, filename, format, options)
  % export as a X3D/XHTML file
  % option can contain: axes (display axes), auto to rescale axes in range [0,1]
  % for aspect ratio 1.
  
  titl = char(a);
  titl(titl=='<')='[';
  titl(titl=='>')=']';
  desc = evalc('disp(a)');
  desc(desc=='<')='[';
  desc(desc=='>')=']';
  a = reducevolume(a);
  a = iData_private_cleannaninf(a);
  x1=min(a); x2=max(a);
  if ndims(a) <= 2
    f=figure('visible','off');
    if ~isempty(strfind(options, 'auto')) || ~isempty(strfind(options, 'tight'))
      setaxis(a,1,linspace(x1,x2,size(a,1)));
      setaxis(a,2,linspace(x1,x2,size(a,2)));
    end
    if ischar(options)
      h = plot(a,options); % make sure file is not too big
    else
      h = plot(a);
    end
    figure2xhtml(filename, f, struct('interactive',true, ...
      'output', strtok(format),'title',titl,'Description',desc));
    close(f);
  else
    [x, xlab] = getaxis(a,2); x=double(x);
    [y, ylab] = getaxis(a,1); y=double(y);
    [z, zlab] = getaxis(a,3); z=double(z);
    [c, clab] = getaxis(a,0); c=double(c);
    ax = 0; iso = [];
    if ischar(options), 
      iso = str2double(strtok(options));
      if strfind(options, 'axes'), ax=1; end
      if ~isempty(strfind(options, 'auto')) || ~isempty(strfind(options, 'tight'))
        y=linspace(x1,x2,size(a,1));
        x=linspace(x1,x2,size(a,2));
        z=linspace(x1,x2,size(a,3));
      end
    end
    if (isempty(iso) || ~isfinite(iso))  
      if isnumeric(options) && isscalar(options)
        iso = options;
      else
        iso = (min(c(:))+max(c(:)))/2;
      end
    end
    fv=isosurface(x,y,z,c,iso);
    x3mesh(fv.faces, fv.vertices, ...
      'format', strtok(format), ...
      'name', filename, 'subheading', desc, 'rotation', 0, 'axes', ax);
    % Create the supporting X3DOM  files
%    folder = fileparts(filename);
%    load(fullfile(fileparts(which('figure2xhtml')), 'functions', 'x3dom.mat'))
%    folder=fullfile(folder, 'x3dom');
%    if(~isdir(folder)), mkdir(folder); end
%    for i=1:length(data)
%        fid = fopen(fullfile(folder, data(i).filename), 'w','ieee-le');
%        fwrite(fid,data(i).filedata,'uint8');
%        fclose(fid);
%    end
  end
