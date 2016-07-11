functions filename = iData_private_saveas_x3d(a, filename)
  titl = char(a);
  titl(titl=='<')='[';
  titl(titl=='>')=']';
  desc = evalc('disp(a)');
  desc(desc=='<')='[';
  desc(desc=='>')=']';
  if ndims(a) <= 2
    f=figure('visible','off');
    if ischar(options)
      h = plot(reducevolume(a),options); % make sure file is not too big
    else
      h = plot(reducevolume(a));
    end
    figure2xhtml(filename, f, struct('interactive',true, ...
      'output', format,'title',titl,'Description',desc));
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
      if strfind(options, 'auto')
        x=linspace(0,1,size(a,1));
        y=linspace(0,1,size(a,2));
        z=linspace(0,1,size(a,3));
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
    x3mesh(fv.faces, fv.vertices, 'name', filename, 'subheading', desc, 'rotation', 0, 'axes', ax);
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
