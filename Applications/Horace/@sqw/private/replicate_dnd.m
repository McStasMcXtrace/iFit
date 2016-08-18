% Make a higher dimensional dataset from a lower dimensional dataset by
% replicating the data along the extra dimensions of a reference dataset.
%
% Syntax:
%   >> dout = replicate (din, dref)
%
% Input:
% ------
%   din     dnd structure.
%
%   dref    Reference dnd structure to use as template for expanding the 
%           input straucture.
%           - The plot axes of din must also be plot axes of dref, and the number
%           of points along these common axes must be the same, although the
%           numerical values of the coordinates need not be the same.
%           - The data is expanded along the plot axes of dref that are 
%           integration axes of din. 
%           - The annotations etc. are taken from the reference dataset.
%
% Output:
% -------
%   dout    Output dnd structure.
%
%