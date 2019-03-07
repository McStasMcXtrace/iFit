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
% other formats are defined by each 'read_*' function, and obtained with 
%   read_FMT('identify')
%
% formats should be sorted from the most specific to the most general.
% Formats will be sorted with text based ones with patterns first, then text based,
% then other formats (e.g. binary).
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
    
    chalkriver.name     ='ChalkRiver CNBC';
    chalkriver.patterns ={'Run ','Seq ','Rec ','Mode ','Temp:','File '};
    chalkriver.options  ='--fast --binary  --headers --comment=NULL --silent --section=Run --metadata=File --metadata=Sig';
    chalkriver.method   ='read_anytext';
    chalkriver.postprocess='load_chalkriver';
    
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
    
% definition of configuration ==================================================
    config.loaders =  { ...
       chalkriver, ...
       qd_vms, ISIS_SQW, agilent_ms, thermo_ms, ...
    };
	       
	  config.UseSystemDialogs = 'yes'; % no: use uigetfiles, else defaults to 'uigetfile'
	  config.FileName         = [ mfilename ' (default configuration from ' which(mfilename) ')' ];
	  
