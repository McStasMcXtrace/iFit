function textart(data,filename)

  % ramp=['`''..-,:"_!~/\;*|^<+7r?v=iJLlYc)T{(}tsIVCxF325]1[uU4nzAXjfoZSyPweaKEHkGOh0M$N9#dq6RmDW%bpQ8Bg@&'];
  ramp = '$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\|()1{}[]?-_+~<>i!lI;:,"^`''. ';
  ramp = ramp(end:-1:1);
  
  % If 3-D image (color) convert to black/white
  if ndims(data) == 3
    data = mean(data,3);
  end
  
  % rescale the image so that its values go from 1 to numel(ramp)
  
  data = data-min(data(:));
  data = data/max(data(:));
  data = floor(data*(numel(ramp)-1))+1;
  
  % now replace the image by its ASCII char from the ramp
  data = ramp(data);

  % write the file
  fid = fopen(filename, 'w');
  for ii=1:size(data,1)
      fprintf(fid, '%s\n', data(ii,:));
  end
  fclose(fid);

