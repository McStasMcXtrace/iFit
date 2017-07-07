function filename = iData_private_saveas_avi(a, filename,options)

% first reduce the object to 2D/3D
if ndims(a) <= 1
  filename = [];
  return; 
elseif ndims(a) > 3
  % reduce dimensions: look for the axes which have the largest extension
  extent = 1:ndims(a);
  for index=1:ndims(a)
    this_axis = getaxis(a, index);
    extent(index) = max(this_axis(:)) - min(this_axis(:));
  end
  % now get the largest extents
  [extent,extent_index] = sort(extent);
  extent_index = extent_index(end:-1:1);  % descending order
  extent       = extent(end:-1:1);
  % we use extend_index(1:3)
  % we keep axes [1:2 ndims(a)] (this allows e.g. to plot S(q,w))
  if isvector(a) > 1 % event data set
    to_remove = [];
    for index=2:ndims(a)
      if isempty([ strfind(method,'whole') strfind(method,'full') ])
        if index <= 4 && extent(index) > extent(1)*1e-2, continue; end
      elseif index <= 3, continue; end
      to_remove = [ to_remove extent_index(index) ];
    end
    a = rmaxis(a,to_remove);
  else
    sz = size(a); sz(extent_index(4:end)) = 1;
    iData_private_warning(mfilename, [ 'Reducing ' num2str(ndims(a)) '-th dimensional data ' a.Tag ' "' a.Title '" to 3D with a=resize(a, ' mat2str(sz) ')' ]);
    a = squeeze(resize(a, sz));
  end
end

if exist('avifile') 
  fig=figure('Visible','off');
  plot(a, options);
  aviobj = avifile(filename);
  for ang=0:5:360
    view([ ang 30])
    aviobj = addframe(aviobj,getframe(fig));
  end
  close(fig)
  aviobj = close(aviobj);
elseif exist('VideoWriter')
  fig=figure('Visible','off');
  plot(a, options);
  aviobj = VideoWriter(filename);
  open(aviobj);
  for ang=0:5:360
    view([ ang 30])
    writeVideo(aviobj,getframe(fig))
  end
  close(aviobj)
  close(fig)
else filename = [];
end

