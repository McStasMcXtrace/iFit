function compiled = muphocor_compile_binary(compile)
  % compile muphocor as binary when does not exist yet
  
  compiled = ''; 
  if nargin == 0, compile = ''; end
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else           precmd=''; end
  
  if ispc, ext='.exe'; else ext=''; end
  this_path = fullfile(fileparts(which(mfilename)));
  
  % try in order: global(system), local, local_arch
  for try_target={ ...
          [ 'muphocor' ext ], 'muphocor'
          fullfile(this_path, [ 'muphocor_' computer('arch') ext ]), ...
          fullfile(this_path, [ 'muphocor' ext ]) }
      
    [status, result] = system([ precmd try_target{1} ]); % run from Matlab

    if status == 1 && nargin == 0
        % the executable is already there. No need to make it .
        compiled = try_target{1};
        return
    end
  end
  
  % when we get there, compile muphocor_arch, not existing yet
  target = fullfile(this_path, [ 'muphocor_' computer('arch') ext ]);

  % search for a FORTRAN compiler
  fc = '';
  for try_fc={getenv('FC'),'gfortran','fc','ifc'}
    if ~isempty(try_fc{1})
      [status, result] = system([ precmd try_fc{1} ]);
      if status == 4 || ~isempty(strfind(result,'no input file'))
        fc = try_fc{1};
        break;
      end
    end
  end
  if isempty(fc)
    if ~ispc
      disp([ mfilename ': ERROR: FORTRAN compiler is not available from PATH:' ])
      disp(getenv('PATH'))
      disp([ mfilename ': You may have to extend the PATH with e.g.' ])
      disp('setenv(''PATH'', [getenv(''PATH'') '':/usr/local/bin'' '':/usr/bin'' '':/usr/share/bin'' ]);');
    end
    error('%s: Can''t find a valid FORTRAN compiler. Install any of: gfortran, ifc or define it as env("FC")\n', ...
    mfilename);
  else
    try
      FCFLAGS = [ '-g -fcheck=all -fbounds-check -fdollar-ok -finit-local-zero ' ...
        '-static -std=legacy -fno-automatic -fdefault-real-8 -fno-align-commons -static-libgfortran' ];
      FCFLAGS = textscan(FCFLAGS, '%s', 'Delimiter',' ');
      FCFLAGS = FCFLAGS{1};
      fprintf(1, '%s: compiling MUPHOCOR binary (using %s)...\n', mfilename, fc);
      cmd={fc, '-O2','-o', target, ...
        fullfile(this_path,'muphocor.f90'),FCFLAGS{:}};
      cmd = sprintf('%s ',cmd{:});
      disp(cmd)
      [status, result] = system([ precmd cmd ]);
      if status == 0
        compiled = target;
      else
        disp(result);
      end
    end
  end

  if isempty(compiled) && ~isempty(compile)
    error('%s: Can''t compile muphocor.f90 binary\n       in %s\n', ...
        mfilename, fullfile(this_path));
  end
  

