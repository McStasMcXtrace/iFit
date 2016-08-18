% Load data from ASCII Tobyfit .par file and returns data in the form, requested by horace
%Usage:
%>> det = get_par(sqw,filename)
%>> det = get_par(sqw,filename,'-array')
%  if varargin ('-array' switch is) present, do not convert into detector structure but return
%  initial array, which is 6xNDet array with NDet equal to number of detectors and the column 
%  meaning correspond to the 
%%   Overloaded methods:
%      sqw/get_par
%      rundata/get_par
%      sqw/get_par
%