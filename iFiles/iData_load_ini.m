function loaders = iData_load_ini
% formats = iData_load_ini
%
% User definitions of specific import formats to be used by iData/load
%
% Each format is specified as a structure with the following fields
%   method:   function name to use, called as method(filename, options...)
%   patterns: list of strings to search in data file. If all found, then method
%             is qualified
%   name:     name of the method/format
%   options:  additional options to pass to the method.
%             If given as a string they are catenated with file name
%             If given as a cell, they are given to the method as additional arguments
%   postprocess: function to call after file import, to assign aliases, ...
%             called as iData=postprocess(iData)
%
% formats should be sorted from the most specific to the most general.
% Formats will be tried one after the other, in the given order.
% System wide loaders are tested after user definitions.
%
% See also: iLoad, save, iData/saveas

    format1.name       ='ILL Data (normal integers)';
    format1.patterns   ={'RRRR','AAAA','FFFF','SSSS','IIII'};
    format1.options    ='--headers --fortran --catenate --fast --binary --makerows=IIII --makerows=FFFF';
    format1.method     ='looktxt';
    format1.postprocess='';
    
    format2.name       ='ILL Data (large integers)';
    format2.patterns   ={'RRRR','AAAA','FFFF','SSSS','JJJJ'};
    format2.options    ='--headers --fortran --catenate --fast --binary --makerows=JJJJ --makerows=FFFF';
    format2.method     ='looktxt';
    format2.postprocess='';
    
    format3.name       ='ILL Data (floats only)';
    format3.patterns   ={'RRRR','AAAA','FFFF','SSSS'};
    format3.options    ='--headers --fortran --catenate --fast --binary --makerows=FFFF';
    format3.method     ='looktxt';
    format3.postprocess='';
    
    format4.name       ='ILL Data (general)';
    format4.patterns   ={'SSSS'};
    format4.options    ='--headers --fortran --catenate --fast --binary --makerows=FFFF --makerows=JJJJ --makerows=IIII';
    format4.method     ='looktxt';
    format4.postprocess='';
    
    format5.name       ='ILL TAS Data (polarized)';
    format5.patterns   ={'POSQE:','PARAM:','DATA_:','LOCAL:','USER_:','PAL'};
    format5.options    ='--headers --section=PARAM --section=VARIA --section=ZEROS --section=DATA -section=POLAN --metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE --metadata=INSTR --metadata=COMND';
    format5.method     ='looktxt';
    format5.postprocess='load_ill_tas'; % load_ill_tas
    
    format6.name       ='ILL TAS Data';
    format6.patterns   ={'POSQE:','PARAM:','DATA_:','LOCAL:','USER_:'};
    format6.options    ='--headers --section=PARAM --section=VARIA --section=ZEROS --section=DATA --metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE --metadata=INSTR --metadata=COMND';
    format6.method     ='looktxt';
    format6.postprocess='load_ill_tas'; % load_ill_tas
    
    format7.name       ='SPEC';
    format7.patterns   ={'#F','#D','#S'};
    format7.options    ='--headers --metadata="#S " --comment= ';
    format7.method     ='looktxt';
    format7.postprocess='';
    
    format8.name       ='McStas Scan output';
    format8.patterns   ={'# Numpoints:','# variables:','# title: Scan of'};
    format8.options    =['--headers --comment= --metadata=variables ' ...
                         '--metadata=xlabel --metadata=ylabel'];
    format8.method     ='looktxt';
    format8.postprocess='load_mcstas_scan';
    
    format9.name       ='McStas 2D monitor';
    format9.patterns   ={'Format: McStas with text headers file.','# type: array_2d'};
    format9.options    = ['--headers --comment= --metadata=variables ' ...
		    '--metadata=Errors --metadata=Events --metadata=xlabel ' ...
		    '--metadata=ylabel --metadata=zlabel --metadata=xylimits'];
    format9.method     ='looktxt';
    format9.postprocess='load_mcstas_2d';
    
    format10.name       ='McStas 1D monitor';
    format10.patterns   ={'Format: McStas with text headers file.','# type: array_1d'};
    format10.options    ='--headers --comment= --metadata=variables --metadata=xlabel --metadata=ylabel';
    format10.method     ='looktxt';
    format10.postprocess='load_mcstas_1d';
    
    format11.name       ='McStas sim file';
    format11.patterns   ={'begin simulation','  Format: McStas'};
    format11.options    ='--headers --comment=';
    format11.method     ='looktxt';
    format11.postprocess='load_mcstas_sim';
    
    format12.name       ='ILL TAS Data (light)';
    format12.patterns   ={'POSQE:','PARAM:','DATA_:','LOCAL:','USER_:'};
    format12.options    ='--headers --section=PARAM --section=DATA --metadata=DATA';
    format12.method     ='looktxt';
    format12.postprocess='';
    
    loaders= { format1, format2, format3, format4, format5, format6, ...
	       format7, format8, format9, format10, format11};
    
