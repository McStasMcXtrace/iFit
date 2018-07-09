function add2imformat(description, ext, isa, info, read, write)
% add2imformat: we add/replace entries to the imformats registry

  formats = imformats;
  
  newFormat.ext   = ext;
  newFormat.isa   = isa;  % @(f)not(isempty(imfinfo(self, f)));
  newFormat.info  = info; % @(f)imfinfo(self, f);
  newFormat.read  = read; % @(f)imread(self, f);
  newFormat.write = write;
  newFormat.alpha = 0;
  newFormat.description = description;
  
  index   = find(strcmp({ formats.description },description), 1);
  if isempty(index) % not registered yet
    imformats('add', newFormat);
  else % replace existing entry
    formats(index) = newFormat;
    imformats(formats);
  end
