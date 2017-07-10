function filename = ResLibCal_GenerateReport(filename)
    if nargin == 0, filename = []; end
    %  look for not-opened views
    toclose=[];
    for dim=1:3
      if isempty(findobj(0, 'Tag',[ 'ResLibCal_View' num2str(dim)]))
        toclose=[ toclose dim ];
      end
    end
    %  list parameters
    %  show view2, view3 and geometry
    [out,f1] = ResLibCal_ViewResolution('',1);
    [out,f2] = ResLibCal_ViewResolution('',2);
    [out,f3] = ResLibCal_ViewResolution('',3);
    views = [f1 f2 f3];
    out = ResLibCal_UpdateViews(out, 'force');
    % dump figures
    if isempty(filename)
      tmpd = tempname;
      if ~isdir(tmpd), mkdir(tmpd); end
      filename = fullfile(tmpd, 'rescal.html');
    else
      [tmpd,f,e] = fileparts(filename);
    end
    print(f1, fullfile(tmpd,'geometry'),'-dpng');
    print(f2, fullfile(tmpd,'resolution2'),'-dpng');
    print(f3, fullfile(tmpd,'resolution3'),'-dpng');
    % close the views that were not there
    close(views(toclose));
    % make up a title
    titl = [ 'ResLibCal ' datestr(now) ];
    if isfield(out, 'EXP') && isfield(out.EXP, 'method')
      titl = [ titl ' ' out.EXP.method ];
    end
    if isfield(out, 'resolution') && iscell(out.resolution) && numel(out.resolution) > 1
      titl = [ titl '. Scan of ' num2str(numel(out.resolution)) ' steps' ];
    end
    % create the document
    
    fid = fopen(filename, 'w+');  % create or append to file
    
    fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n');
    fprintf(fid, '<html>\n<head>\n<title>%s</title>\n</head>\n', ...
        titl);
    fprintf(fid, '<body>\n');
    fprintf(fid, '<h2>%s</h2>\n', titl);
    
    % add the Table of Contents
    fprintf(fid, '<hr>Table of contents:<br><ul>\n');
    fprintf(fid, '<li><a href="#Instrument_parameters">Instrument parameters</a></li>\n');
    fprintf(fid, '<li><a href="#Resolution">Resolution</a></li>\n');
    fprintf(fid, '<li><a href="#Views">Views (2D/3D/Geometry)</a></li>\n');
    fprintf(fid, '</ul><hr>\n');
    
    fprintf(fid, 'This document contains the configuration and resolution calculation for a Triple-Axis neutron spectrometer, generated with <a href="http://ifit.mccode.org/Applications/ResLibCal/doc/ResLibCal.html">ResLibCal</a>. Refer to this link in order to get more information on the calculation methods and usage.\n');
    
    % list parameters
		rlu = get(ResLibCal_fig('View_ResolutionRLU'), 'Checked');    % [a* b*  c* ]
		spec= get(ResLibCal_fig('View_ResolutionSPEC'),'Checked');    % [Ql Qt  Qv ]
		abc = get(ResLibCal_fig('View_ResolutionABC'), 'Checked');    % [A  B   C  ]
		lat = get(ResLibCal_fig('View_ResolutionLattice'), 'Checked');% [a* b'* c'*]
		modev='abc'; % default
		if     strcmp(rlu, 'on') modev='rlu'; 
		elseif strcmp(spec,'on') modev='spec'; 
		elseif strcmp(abc, 'on') modev='abc';
		elseif strcmp(lat, 'on') modev='lattice'; end
		[res, inst] = ResLibCal_FormatString(out, modev);
		fprintf(fid, '<h3><a name = "Instrument_parameters">Instrument parameters</h3>\n<pre>\n');
		fprintf(fid, '%s\n', inst{:});
		
		fprintf(fid,'</pre><br>\n');
		fprintf(fid, '<h3><a name = "Resolution">Resolution</h3>\n<pre>\n');
		for l=res'
		  b= l{1};
		  c = cellstr(l{1});
		  fprintf(fid, '%s\n', c{:}); % so that matrices appear OK
		end
    fprintf(fid,'</pre><br>\n');
    
    % add images
    fprintf(fid, '<h3><a name="Views">Resolution and geometry views</h3>\n');
    for name={'geometry','resolution2','resolution3'}
      fprintf(fid, '<div style="text-align: center;"><a href="%s"><img src="%s" align="middle"></a><br>\n<i>View: %s</i><br></div>\n', ...
      [ name{1} '.png' ], ...
      [ name{1} '.png' ], name{1});
    end
    
    % display a 'footer' below the object description
    fprintf(fid,[ '<br><hr><b>' datestr(now) '</b> - ' version(iData) '<br>\n' ]);
    
    fprintf(fid,[ '<a href="http://ifit.mccode.org">Powered by iFit ' ...
      '<img src="http://ifit.mccode.org/images/iFit-logo.png" width=35 height=32></a> \n' ...
      '<a href="http://www.ill.eu">(c) ILL ' ...
      '<img title="ILL, Grenoble, France www.ill.eu" src="http://ifit.mccode.org/images/ILL-web-jpeg.jpg" alt="ILL, Grenoble, France www.ill.eu" style="width: 33px; height: 32px;"></a><hr>\n' ]);
    fprintf(fid, '</body>');
    web(filename);
    if isdeployed
      disp([ mfilename ': Generated report as ' filename ]);
    else
      disp([ mfilename ': Generated report as <a href="' filename '">' filename '</a>' ]);
    end
