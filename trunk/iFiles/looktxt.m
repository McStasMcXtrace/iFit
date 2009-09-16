% [structs] = looktxt('filename [options]') Import text data
% Usage: looktxt [options] file1 file2 ...
% Action: Search and export numerics in a text/ascii file.
%    This program analyses files looking for numeric parts
%    Each identified numeric field is named and exported
%    into a structure with fields
%      ROOT.SECTION.FIELD = VALUE
%    In order to sort your data, you may specify as many --section
%    and --metadata options as necessary
% Return: a single structure or a cell of structures
%    If the structure can not be evaluated, the raw Matlab script is returned.
% Example: looktxt -c -s PARAM -s DATA filename
%
% Useful options when used from Matlab:
% --binary   or -b    Stores numerical matrices into an additional binary
%                     float file, which makes further import much faster.
% --catenate or -c    Catenates similar numerical fields
% --force    or -F    Overwrites existing files
% --fortran --wrapped Catenates single Fortran-style output lines with
%                     previous matrices
% --headers  or -H    Extracts headers for each numerical field
% --section=SEC       Classifies fields into section matching word SEC
%       -s SEC
% --metadata=META     Extracts lines containing word META as user meta data
%         -m META
% --fast              Uses a faster reading method, requiring numerics
%                     to be separated by \n\t\r\f\v and spaces only
% --makerows=NAME     All fields matching NAME are transformed into row vectors
%
% Usual options are: --fast --fortran --binary --force --catenate --comment=NULL
% List of all options can be obtained using: looktxt --help
%
% looktxt  version 1.0.8 (16 Sept 2009) by Farhi E. [farhi@ill.fr]

% if we come here it means looktxt was not made a MeX file
mex -O -output looktxt texmex.c

