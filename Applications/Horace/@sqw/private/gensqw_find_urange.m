% Get limits of data for grid on which to store sqw data
% ---------------------------------------------------------
% Use the fact that the lowest and highest energy transfer bin centres define the maximum extent in
% momentum and energy. We calculate using the full detector table i.e. do not account for masked detectors
% but this is reasonable so long as not many detecotrs are masked. 
% (*** In more systematic cases e.g. spe file is for MARI, and impose a mask file that leaves only the
%  low angle detectors, then the calculation will be off. Will bw able to rectify this once use
%  libisis run file structure, when can enquire of masked detectors from the IXTrunfile object)
%