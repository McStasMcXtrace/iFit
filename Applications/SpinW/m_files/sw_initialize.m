function sw_initialize()
% initializes sw library and cleans symmetry.dat file from user entries
%
% SW_INITIALIZE()
%
% See also SW.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

%cd(sw_rootdir);
%m2html('mfiles','m_files','htmldir','doc')
%builddocsearchdb([sw_rootdir 'doc']);

% clear symmetry.dat file
symPath = [sw_rootdir 'dat_files' filesep 'symmetry.dat'];
% Count the number of lines
fid = fopen(symPath);
fidNew = fopen([symPath '0'],'w');


if fid == -1
    error('spinw:sw_gensym:FileNotFound',['Symmetry definition file not found: '...
        regexprep(symPath,'\' , '\\\') '!']);
end

nLines = 1;
while (~feof(fid)) && (nLines <= 230)
    line = fgetl(fid);
    if nLines < 230
        fprintf(fidNew,[line '\n']);
    else
        fprintf(fidNew,line);
    end
    nLines = nLines + 1;
end
fclose(fid);
fclose(fidNew);

% copy file
delete(symPath);
copyfile([symPath '0'],symPath);
delete([symPath '0']);

end