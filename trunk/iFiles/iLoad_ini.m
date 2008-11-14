function loaders = iLoad_ini
% formats = iLoad_ini
%
% User definitions of specific import formats to be used by iLoad
%
% Each format is specified as a structure with the following fields
%   method:   function name to use, called as method(filename, options...)
%   extension:a single or a cellstr of extensions associated with the method
%   patterns: list of strings to search in data file. If all found, then method
%             is qualified
%   name:     name of the method/format
%   options:  additional options to pass to the method.
%             If given as a string they are catenated with file name
%             If given as a cell, they are given to the method as additional arguments
%   postprocess: function called from iData/load after file import, to assign aliases, ...
%             called as iData=postprocess(iData)
%
% formats should be sorted from the most specific to the most general.
% Formats will be tried one after the other, in the given order.
% System wide loaders are tested after user definitions.
%
% See also: iLoad, save, iData/saveas

    format1.name       ='ILL Data (normal integers)';
    format1.patterns   ={'RRRR','AAAA','FFFF','IIII'};
    format1.options    ='--headers --fortran --catenate --fast --binary --makerows=IIII --makerows=FFFF --silent';
    format1.method     ='looktxt';
    
    format2.name       ='ILL Data (large integers)';
    format2.patterns   ={'RRRR','AAAA','FFFF','JJJJ'};
    format2.options    ='--headers --fortran --catenate --fast --binary --makerows=JJJJ --makerows=FFFF --silent';
    format2.method     ='looktxt';
    
    format3.name       ='ILL Data (floats only)';
    format3.patterns   ={'RRRR','AAAA','FFFF'};
    format3.options    ='--headers --fortran --catenate --fast --binary --makerows=FFFF --silent';
    format3.method     ='looktxt';
    
    format4.name       ='ILL Data (general)';
    format4.patterns   ={'SSSS'};
    format4.options    ='--headers --fortran --catenate --fast --binary --makerows=FFFF --makerows=JJJJ --makerows=IIII --silent';
    format4.method     ='looktxt';
    
    format5.name       ='ILL TAS Data (polarized)';
    format5.patterns   ={'POSQE:','PARAM:','DATA_:','LOCAL:','USER_:','PAL'};
    format5.options    =['--fast --binary --headers ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=POLAN --section=STEPS ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE'];
    format5.method     ='looktxt';
    format5.postprocess='load_ill_tas'; % load_ill_tas
    format5.extension  ='scn';
    
    format6.name       ='ILL TAS Data';
    format6.patterns   ={'POSQE:','PARAM:','DATA_:','LOCAL:','USER_:'};
    format6.options    =['--fast --binary --headers ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=STEPS ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE'];
    format6.method     ='looktxt';
    format6.postprocess='load_ill_tas'; % load_ill_tas
    format6.extension  ='scn';
    
    format7.name       ='SPEC';
    format7.patterns   ={'#F','#D','#S'};
    format7.options    ='--fast --binary --headers --metadata="#S " --comment= ';
    format7.method     ='looktxt';
    format7.extension  ='spc';
    
    format8.name       ='McStas Scan output';
    format8.patterns   ={'# Numpoints:','# variables:','# title: Scan of'};
    format8.options    =['--fast --binary --headers --comment= --metadata=variables ' ...
                         '--metadata=xlabel --metadata=ylabel'];
    format8.method     ='looktxt';
    format8.postprocess='load_mcstas_scan';
    
    format9.name       ='McStas 2D monitor';
    format9.patterns   ={'Format: McStas with text headers file.','# type: array_2d'};
    format9.options    = ['--fast --binary --headers --comment= --metadata=variables ' ...
		    '--metadata=Errors --metadata=Events --metadata=xlabel ' ...
		    '--metadata=ylabel --metadata=zlabel --metadata=xylimits'];
    format9.method     ='looktxt';
    format9.postprocess='load_mcstas_2d';
    
    format10.name       ='McStas 1D monitor';
    format10.patterns   ={'Format: McStas with text headers file.','# type: array_1d'};
    format10.options    ='--fast --binary --headers --comment= --metadata=variables --metadata=xlabel --metadata=ylabel';
    format10.method     ='looktxt';
    format10.postprocess='load_mcstas_1d';
    
    format11.name       ='McStas sim file';
    format11.extension  ='sim';
    format11.patterns   ={'begin simulation','  Format: McStas'};
    format11.options    ='--fast --binary --headers --comment=';
    format11.method     ='looktxt';
    format11.postprocess='load_mcstas_sim';
    
    format13.name       ='McStas Sqw table';
    format13.patterns   ={'Sqw data file for Isotropic_Sqw'};
    format13.options    ='--fast --binary --headers';
    format13.method     ='looktxt';
    format13.postprocess='load_mcstas_sqw';
    format13.extension  ='sqw';
    
    format14.name       ='INX tof data';
    format14.options    ='--headers --fortran --catenate --fast --binary';
    format14.method     ='looktxt';
    format14.postprocess='load_ill_inx';
    format14.patterns   ={'INX'};
    format14.extension  ='inx';
    
    loaders= { format1, format2, format3, format4, format5, format6, ...
	       format7, format8, format9, format10, format11, format13, format14 };
    
