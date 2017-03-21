function filename = iData_private_saveas_stl(a, filename, format)
  if ndims(a) == 1    iData_private_warning(mfilename,[ 'Object ' inputname(1) ' ' a.Tag ' 1D does not seem to be exportatble as a ' format ' file. Ignoring.' ]);
    return
  else
    a = iData_private_cleannaninf(a);
    if ndims(a) == 2
      [x, xlab] = getaxis(a,2); x=double(x);
      [y, ylab] = getaxis(a,1); y=double(y);
      [z, zlab] = getaxis(a,0); z=double(z);
    elseif ndims(a) >= 3
      [x, xlab] = getaxis(a,2); x=double(x);
      [y, ylab] = getaxis(a,1); y=double(y);
      [z, zlab] = getaxis(a,3); z=double(z);
    else
      disp([ mfilename ': Export into ' format ' is only possible for 2D+ objects, not for ' num2str(ndims(a)) 'D. Ignoring.' ])
      return
    end
    if any(strcmp(format, {'stl','stlb'}))
      mode = 'binary';
    else
      mode = 'ascii';
    end
    
    % get Title
    T   = a.Title; if ~ischar(T), T=char(T); end
    if ~isvector(T), T=transpose(T); T=T(:)'; end
    T   = regexprep(T,'\s+',' '); % remove duplicated spaces
    if length(T) > 69, T=[ T(1:60) '...' T((end-8):end) ]; end
    % get the faces and vertices
    if isvector(x) && isvector(y)
      [x,y] = meshgrid(x,y);
    end
    v = [x(:) y(:) z(:)];
    f = delaunay(x(:),y(:));
    %f=delaunayn(v,{'Qx','Qv','Tv', 'Qt','Qbb','Qc','Qz'});
    [v, indexm, indexn] =  unique(v, 'rows');
    f = indexn(f);
    if strncmp(format,'stl',3)  % STL format
      stlwrite(filename, f,v, 'Mode', mode, 'Title', T);
    else                        % OFF and PLY formats
      [fid, message]=fopen(filename,'w+');
      if fid == -1
        iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag 'in format ' format ]);
        disp(message)
        return
      end
      % write header
      if strcmp(format,'ply')   % OFF format
        fprintf(fid,'ply\nformat ascii 1.0\ncomment This is a PLY file format for %s\nelement vertex %d\n', ...
          T, size(v,1));
        fprintf(fid,'property float x\nproperty float y\nproperty float z\nelement face %d\n', size(f,1));
        fprintf(fid,'property list uchar int vertex_indices\nend_header\n');
      else                      % PLY
        fprintf(fid,'OFF\n%i %i 0\n',...
          size(v,1), size(f,1));
      end
      % write data
      str = num2str(v,5); str(:,end+1) = sprintf('\n');
      str = str'; str = str(:)';
      fprintf(fid,'%s', str);
      f = [ (size(f,2)+1)*ones(size(f,1),1) f ];
      str = num2str(f-1); str(:,end+1) = sprintf('\n');
      str = str'; str = str(:)';
      fprintf(fid,'%s', str);
      if strcmp(format,'off')
        fprintf(fid,'# This is an Object File Format (geomview) for %s\n', T);
       end
      fclose(fid);
    end
  end
