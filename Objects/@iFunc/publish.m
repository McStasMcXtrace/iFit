function a = publish(a, filename)
% b = publish(s) : export the Model object as an HTML document.
%
%   @iFunc/publish function to export an object into a readable document.
%
%   publish(a, file)  export to given file name
%   publish(a, dir)   export to given directory
%
% input:  s: object or array (iFunc)
%         filename: a file name or directory where to export. When not given
%           the exportation takes place in the local directory.
% output: b: generated filename
% ex:     b=publish(gauss); webbrowser(b, 'system')
%
% Version: $Date$
% See also iFunc, iFunc/save

if nargin < 2, filename = ''; end

if ndims(a) == 4 && strncmp(a.Name, 'Sqw_Phonon', numel('Sqw_Phonon'))
  a = sqw_phonons(a, 'report');
  if isfield(a.UserData,'dir') && isdir(a.UserData.dir) ...
    && ~isempty(dir(fullfile(a.UserData.dir, 'sqw_phonons.html')))
    a = fullfile(a.UserData.dir, 'sqw_phonons.html');
  end
else
  a = saveas(a, filename, 'html');
end

if nargout == 0
  web(a)
end

