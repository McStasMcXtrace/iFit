% Get data from ascii file with column data qx-qy-qz-eps-signal-error
%
%   >> data = get_ascii_column_data (datafile)
%
% Input:
% ------
%   datafile    Full file name of ascii data file
%               Format is one of the following column arrangements:
%                   qx'  qy'  qz'  S
%                   qx'  qy'  qz'  S  ERR
%                   qx'  qy'  qz'  eps  S  ERR
%
%               Here qz' is the component of momentum along ki (Ang^-1)
%                    qy' is component vertically upwards (Ang^-1)
%                    qx' defines a hight-hand coordinate frame with qy' and qz'
%                    S   signal
%                    ERR standard deviation
%
% Output:
% -------
%   data        Data structure with following fields:
%                   data.filename   Name of file excluding path
%                   data.filepath   Path to file including terminating file separator
%                   data.qspec      [4 x n] array of qx,qy,qz,eps of all the data points
%                                  where now the component are in spectrometer coordinates
%                                 (qx||ki, qz up; qx,qy,qz orthonormal and units Ang^-1)
%                   data.S          [1 x n] array of signal values
%                   data.ERR        [1 x n] array of error values (st. dev.)
%                   data.en         Column vector length 2 of min and max eps in the ascii file
%   det         Data structure containing fake detector parameters for unmasked
%              detectors (see get_par for fields)
%