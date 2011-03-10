function config = iLoad_ini
% config = iLoad_ini
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
% These formats can be obtained using [config, configfile]=iLoad('','load config').
% the iLoad_ini configuration file can be saved in the Preference directory
% using [config, configfile] = iLoad(config,'save config').
% A list of all supported formats is shown with iLoad('formats');
%
% See also: iLoad, save, iData/saveas

% definition of formats
    ILL_normal.name       ='ILL Data (normal integers)';
    ILL_normal.patterns   ={'RRRR','AAAA','FFFF','IIII'};
    ILL_normal.options    ='--headers --fortran --catenate --fast --binary --makerows=IIII --makerows=FFFF --silent';
    ILL_normal.method     ='looktxt';
    
    ILL_integers.name       ='ILL Data (large integers)';
    ILL_integers.patterns   ={'RRRR','AAAA','FFFF','JJJJ'};
    ILL_integers.options    ='--headers --fortran --catenate --fast --binary --makerows=JJJJ --makerows=FFFF --silent';
    ILL_integers.method     ='looktxt';
    
    ILL_float.name       ='ILL Data (floats only)';
    ILL_float.patterns   ={'RRRR','AAAA','FFFF'};
    ILL_float.options    ='--headers --fortran --catenate --fast --binary --makerows=FFFF --silent';
    ILL_float.method     ='looktxt';
    
    ILL_general.name       ='ILL Data (general)';
    ILL_general.patterns   ={'SSSS'};
    ILL_general.options    ='--headers --fortran --catenate --fast --binary --makerows=FFFF --makerows=JJJJ --makerows=IIII --silent';
    ILL_general.method     ='looktxt';
    
    ILL_TAS_pol.name       ='ILL TAS Data (polarized)';
    ILL_TAS_pol.patterns   ={'PAL','POSQE:','PARAM:','DATA_:','LOCAL:','USER_:'};
    ILL_TAS_pol.options    =['--fast --binary --headers --silent ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=POLAN --section=STEPS ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE --metadata=MULTI '];
    ILL_TAS_pol.method     ='looktxt';
    ILL_TAS_pol.postprocess='load_ill_tas'; % load_ill_tas
    ILL_TAS_pol.extension  ='scn';
    
    ILL_TAS.name       ='ILL TAS Data';
    ILL_TAS.patterns   ={'POSQE:','PARAM:','DATA_:','LOCAL:','USER_:'};
    ILL_TAS.options    =['--fast --binary --headers --silent ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=STEPS ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE --metadata=DATA ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE --metadata=MULTI '];
    ILL_TAS.method     ='looktxt';
    ILL_TAS.postprocess='load_ill_tas'; % load_ill_tas
    ILL_TAS.extension  ='scn';
    
    spec.name       ='SPEC';
    spec.patterns   ={'#F','#D','#S'};
    spec.options    ='--fast --binary --headers --metadata="#S " --comment=NULL --silent';
    spec.method     ='looktxt';
    spec.extension  ='spc';
    
    mcstas_scan.name       ='McStas Scan DAT output';
    mcstas_scan.patterns   ={'# type: multiarray_1d','# variables:','# title: Scan of'};
    mcstas_scan.options    =['--fast --binary --headers --comment=NULL --metadata=variables --silent ' ...
                         '--metadata=xlabel --metadata=ylabel --metadata=xvars' ];
    mcstas_scan.method     ='looktxt';
    mcstas_scan.postprocess='load_mcstas_scan';
    
    mcstas_2D.name       ='McStas 2D monitor';
    mcstas_2D.patterns   ={'Format: McStas with text headers','# type: array_2d'};
    mcstas_2D.options    = ['--fast --binary --headers --comment=NULL --metadata=variables --silent ' ...
		    '--metadata=Errors --metadata=Events --metadata=xlabel ' ...
		    '--metadata=ylabel --metadata=zlabel --metadata=xylimits --metadata=component' ];
    mcstas_2D.method     ='looktxt';
    mcstas_2D.postprocess='load_mcstas_2d';
    
    mcstas_1D.name       ='McStas 1D monitor';
    mcstas_1D.patterns   ={'Format: McStas with text headers','# type: array_1d'};
    mcstas_1D.options    =['--fast --binary --headers --comment=NULL --silent --metadata=variables ' ...
        '--metadata=xlabel --metadata=ylabel  --metadata=component' ];
    mcstas_1D.method     ='looktxt';
    mcstas_1D.postprocess='load_mcstas_1d';
    
    mcstas_sim.name       ='McStas sim file';
    mcstas_sim.extension  ='sim';
    mcstas_sim.patterns   ={'begin simulation','Format: McStas'};
    mcstas_sim.options    ='--fast --binary --headers  --comment=NULL --silent';
    mcstas_sim.method     ='looktxt';
    mcstas_sim.postprocess='load_mcstas_sim';
    
    mcstas_sqw.name       ='McStas Sqw table';
    mcstas_sqw.patterns   ={'Sqw data file for Isotropic_Sqw'};
    mcstas_sqw.options    ='--fast --binary  --headers --comment=NULL --silent';
    mcstas_sqw.method     ='looktxt';
    mcstas_sqw.postprocess='load_mcstas_sqw';
    mcstas_sqw.extension  ='sqw';
    
    ISIS_spe.name       ='ISIS/SPE tof data';
    ISIS_spe.options    ='--headers --fortran  --catenate --fast --binary --comment=NULL --silent';
    ISIS_spe.method     ='looktxt';
    ISIS_spe.postprocess='load_ill_spe';
    ISIS_spe.patterns   ={'Phi Grid'};
    ISIS_spe.extension  ='spe';
    
    ILL_inx.name       ='INX tof data';
    ILL_inx.options    ='--headers --fortran  --catenate --fast --binary --comment=NULL --silent';
    ILL_inx.method     ='looktxt';
    ILL_inx.postprocess='load_ill_inx';
    ILL_inx.patterns   ={'INX'};
    ILL_inx.extension  ='inx';
    
    ESRF_edf.name       ='EDF ESRF Data Format';
    ESRF_edf.options    ='';
    ESRF_edf.method     ='medfread';
    ESRF_edf.extension  ='edf';
    
% definition of configuration
    config.loaders =  { ILL_normal, ILL_integers, ILL_float, ILL_general, ILL_TAS_pol, ILL_TAS, ...
	       spec, mcstas_scan, mcstas_2D, mcstas_1D, mcstas_sim, mcstas_sqw, ISIS_spe, ILL_inx, ESRF_edf };
	       
	  config.UseSystemDialogs = 'no'; % no: use uigetfiles, else defaults to 'uigetfile'
	  config.FileName         = [ mfilename ' (default configuration from ' which(mfilename) ')' ];
    
