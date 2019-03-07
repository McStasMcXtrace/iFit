function config = iLoad_ini
% config = iLoad_ini User definitions of specific import formats to be used by iLoad
%
% Each format is specified as a structure with the following fields
%   method:   function name to use, called as method(filename, options...)
%   extension:a single or a cellstr of extensions associated with the method
%   patterns: list of strings to search in data file. If all found, then method
%             is qualified. The patterns can be regular expressions.
%             When given as a string, the file is assumed to contain "text" data.
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
% Example: config = iLoad_ini; isstruct(config)
%
% See also: iLoad, save, iData/saveas
% (c) E.Farhi, ILL. License: EUPL.

% definition of formats ========================================================

    ILL_normal.name       ='ILL Data (normal integers)';
    ILL_normal.patterns   ={'RRRR','AAAA','FFFF','IIII'};
    ILL_normal.options    ='--headers --fortran --catenate --fast --binary --makerows=IIII --makerows=FFFF --silent ';
    ILL_normal.method     ='read_anytext';
    
    ILL_integers.name       ='ILL Data (large integers)';
    ILL_integers.patterns   ={'RRRR','AAAA','FFFF','JJJJ'};
    ILL_integers.options    ='--headers --fortran --catenate --fast --binary --makerows=JJJJ --makerows=FFFF --silent ';
    ILL_integers.method     ='read_anytext';
    
    ILL_float.name          ='ILL Data (floats only)';
    ILL_float.patterns      ={'RRRR','AAAA','FFFF'};
    ILL_float.options       ='--headers --fortran --catenate --fast --binary --makerows=FFFF --silent ';
    ILL_float.method        ='read_anytext';
    
    ILL_general.name       ='ILL Data (general)';
    ILL_general.patterns   ={'RRRR','AAAA','SSSS'};
    ILL_general.options    ='--headers --fortran --catenate --fast --binary --makerows=FFFF --makerows=JJJJ --makerows=IIII --silent ';
    ILL_general.method     ='read_anytext';
    
    HZB_FELXX_Flat.name    ='HZB FLEXX FlafCone';
    HZB_FELXX_Flat.patterns={'RRRR','AAAA','VVVV','DATA_:','flat:','POSQE:','PARAM:'};
    HZB_FELXX_Flat.options =[ '--headers --fortran --catenate --fast --binary --silent ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=STEPS --metadata=KFIX ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE --metadata=flat --metadata=PNT'];
    HZB_FELXX_Flat.method  ='read_anytext';
    
    ILL_TAS_pol.name       ='ILL TAS Data (polarized)';
    ILL_TAS_pol.patterns   ={'PAL','POSQE:','PARAM:','DATA_:','USER_:'};
    ILL_TAS_pol.options    =['--fast --binary --headers --silent --fortran=0 --catenate ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=POLAN --section=STEPS ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE --metadata=MULTI --metadata=PNT'];
    ILL_TAS_pol.method     ='read_anytext';
    ILL_TAS_pol.postprocess='load_ill_tas'; % load_ill_tas
    ILL_TAS_pol.extension  ='scn';
    
    ILL_TAS.name       ='ILL TAS Data';
    ILL_TAS.patterns   ={'POSQE:','PARAM:','DATA_:','USER_:'};
    ILL_TAS.options    =['--fast --binary --headers --silent --fortran=0 --catenate ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=STEPS ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE --metadata=DATA ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE --metadata=MULTI --metadata=PNT '];
    ILL_TAS.method      ='read_anytext';
    ILL_TAS.postprocess ='load_ill_tas'; % load_ill_tas
    ILL_TAS.extension   ='scn';
    
    spec.name           ='SPEC';
    spec.patterns       ={'#F','#D','#S'};
    spec.options        ='--fast --binary --headers --metadata=''#S '' --comment=NULL --silent ';
    spec.method         ='read_anytext';
    spec.extension      ={'spc','spec'};
    
    chalkriver.name     ='ChalkRiver CNBC';
    chalkriver.patterns ={'Run ','Seq ','Rec ','Mode ','Temp:','File '};
    chalkriver.options  ='--fast --binary  --headers --comment=NULL --silent --section=Run --metadata=File --metadata=Sig';
    chalkriver.method   ='read_anytext';
    chalkriver.postprocess='load_chalkriver';
    
    CFL.name            ='CFL FullProf crystallography file';
    CFL.patterns        ={'Spgr','Atom'};
    CFL.method          ='read_anytext';
    CFL.options         ='--headers --fortran --catenate --fast --binary --section=Atom --silent --metadata=Spgr --metadata=Cell --metadata=Spgr ';
    CFL.extension       ={'cfl'};
    CFL.postprocess     ='opencfl';
    
    xye.name            ='Simple x/y/e column file';
    xye.extension       ={'xye', 'xy'};
    xye.method          ='read_anytext';
    xye.patterns        ={};
    xye.postprocess     = 'load_xyen';

    OFF_ascii.name      ='OFF 3D ascii';
    OFF_ascii.method    ='read_anytext';
    OFF_ascii.options   ='--fast --binary --headers --comment=NULL --metadata=OFF --silent';
    OFF_ascii.extension ='off';
    OFF_ascii.patterns  ={'\<OFF\>'};
    OFF_ascii.postprocess='openoff';
    
    PLY_ascii.name      ='PLY 3D ascii';
    PLY_ascii.method    ='read_anytext';
    PLY_ascii.options   ='--fast --binary --headers --comment=NULL --silent';
    PLY_ascii.extension ='ply';
    PLY_ascii.patterns  ={'ply','format ascii','element','end_header'};
    PLY_ascii.postprocess='openply';
    
    EZD.name            ='EZD electronic density map';
    EZD.method          ='read_anytext';
    EZD.options         ='--fortran --headers --binary --fast --catenate --comment=NULL --silent';
    EZD.extension       ='ezd';
    EZD.patterns        ={'EZD_MAP','CELL','EXTENT'};
    EZD.postprocess     ='this.Data.MAP = reshape(this.Data.MAP, this.Data.EXTENT);';
    
    qd_vms.name         = 'Quantum Design VMS ppms/mpms';
    qd_vms.extension    = 'dat';
    qd_vms.method       = 'read_anytext';
    qd_vms.options      = '-H -m Headers -m INFO -m Data -m mm --catenate --fast --binary --silent';
    qd_vms.patterns     = {'Quantum Design'};
    
% binary formats ===============================================================
    
    ISIS_SQW.name       = 'ISIS Horace SQW';
    ISIS_SQW.method     = 'sqw';
    ISIS_SQW.extension  = 'sqw';
    ISIS_SQW.postprocess= 'this.Data = struct(this.Data); this = setaxis(this,0);';
    
    agilent_ms.name     = 'Agilent Mass Spectrometry/Chromatography LC/MS GC/MS GC/FID';
    agilent_ms.method   = 'ImportAgilent';  % private
    agilent_ms.extension= {'ch','ms','d'};
    
    thermo_ms.name      = 'Thermo Finnigan Mass Spectrometry/Chromatography';
    thermo_ms.method    = 'ImportThermo'; % private
    thermo_ms.extension = 'raw';
    
% definition of configuration
    config.loaders =  { ILL_normal, ILL_integers, ILL_float, ILL_general, ILL_TAS_pol, ILL_TAS, ...
      spec, chalkriver, HZB_FELXX_Flat, OFF_ascii, PLY_ascii, CFL, EZD, ...
      xye, qd_vms, ISIS_SQW, agilent_ms, thermo_ms, ...
    };
	       
	  config.UseSystemDialogs = 'yes'; % no: use uigetfiles, else defaults to 'uigetfile'
	  config.FileName         = [ mfilename ' (default configuration from ' which(mfilename) ')' ];
	  
