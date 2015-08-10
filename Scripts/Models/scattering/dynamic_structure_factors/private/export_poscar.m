function [ status ] = export_poscar( filename, geometry, option )
%EXPORT_POSCAR Export a geometry struct as a VASP POSCAR file.
%   status = EXPORT_POSCAR(filename,geometry) exports the geometry as a 
%   VASP POSCAR file. 
%
%   EXPORT_POSCAR(filename,geometry,'nosymbols') does not write the elements 
%     line before the coordinates
%
%   See also IMPORT_POSCAR.
%
% VaspLab (BSD) by Max Radin 23 May 2012 (Updated 10 Dec 2012) 
% http://www.mathworks.com/matlabcentral/fileexchange/36836-vasplab

    % write selective dynamics and chemical symbols

    fid = fopen(filename,'w');
    if fid==-1
        error([mfilename ': Error opening ' filename]); 
    end
    if nargin < 3, option=''; end
    
    fprintf(fid,[geometry.comment '\n']);
    fprintf(fid,'1.0\n'); % scale factor
    fprintf(fid, '%19.16f %19.16f %19.16f\n', geometry.lattice'); % lattice vectors
    
    % this block can be unactvated to stay 100% compatible with VASP and PHON
    if ~isempty(geometry.symbols) && (isempty(option) || ~strcmp(option, 'nosymbols'))
        cellfun(@(x) fprintf(fid, '%s ', x), geometry.symbols);
        fprintf(fid, '\n');
    end
    
    fprintf(fid, '%d ', geometry.atomcount); % number of each species
    fprintf(fid, '\nDirect\n');
    fprintf(fid, '%19.16f %19.16f %19.16f\n', geometry.coords');   
    
    fclose(fid);

    status = 0;
    
end

