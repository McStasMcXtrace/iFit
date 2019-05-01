function r = runtest(s, m)
% RUNTEST run a set of tests on object methods

if nargin == 1
  m = methods(s);
end
m=cellstr(m);

% perform all tests
r = runtest_dotest(s, m);

% display results
fprintf(1, '%s: %s %i tests\n', mfilename, class(s), numel(r));
disp('--------------------------------------------------------------------------------')
runtest_report(r);

% ------------------------------------------------------------------------------
function r = runtest_dotest(s, m)
% RUNTEST_DOTEST performs a test for given method 'm'
  r = [];
  fprintf(1,'Running %s (%s) for %i methods\n',...
    mfilename, class(s), numel(m));

  for mindex=1:numel(m)
    % test if the method is functional
    if strcmp(mfilename, m{mindex}), continue; end % not myself
    res.Name      =m{mindex};
    res.Passed    =false;
    res.Failed    =false;
    res.Incomplete=false;
    res.Duration  =0;
    res.Details.nargin =0; 
    res.Details.nargout=0; 
    res.Details.Code   =''; 
    res.Details.Index  =0;
    res.Details.Status ='';
    res.Details.Output ='';
    
    try
      fun = eval([ '@s.' m{mindex} ]); 
      res.Details.nargin = nargin(fun);
      res.Details.nargout= nargout(fun);
      
      h = help([ class(s) '.' m{mindex} ]);
    catch ME
      res.Details.Output = ME;
      h=[];
    end
    if isempty(h)
      res.Details.Status = 'invalid';
      res.Failed         = true;
      res.Incomplete     = true;
      res.Details.Code   = fun;
      r = [ r res ];
      continue
    end
    h = textscan(h, '%s','Delimiter','\n'); h=strtrim(h{1});
    % get the Examples: lines
    ex= strfind(strtrim(h), 'Example:');
    for index=find(~cellfun(@isempty, ex))
      if ~isempty(index) && ~isempty(h{index})
        res.Details.Code = strrep(h{index}, 'Example:','');
        res.Details.Index= index; % line index in 'help'
      else res.Details.Code = ''; end
      if isempty(res.Details.Code)
        res.Details.Status = 'notest';
        res.Passed         = true;
        r = [ r res ];
        break; 
      end
      
      t0 = clock;
      [passed, failed, output] = runtest_sandbox(res.Details.Code);
      res.Passed = passed;
      res.Failed = failed;
      res.Details.Output = output;
      if passed,     res.Details.Status = 'passed';
      elseif failed, res.Details.Status = 'FAILED';
      end
      res.Duration = etime(clock, t0);
      r = [ r res ];
    end
    
  end
  fprintf(1, '\n');
  disp([ 'Done ' mfilename ' (' class(s) ')' ])
% ------------------------------------------------------------------------------
function runtest_report(r)

  % we display details for Incomplete and Failed tests
  for index=1:numel(r)
    res = r(index);
    if res.Failed || res.Incomplete
      % get the error message
      if isa(res.Details.Output,'MException')
        ME = res.Details.Output;
        if ~isempty(ME.stack) && ~strcmp(ME.stack(1).name, 'runtest_sandbox')
          msg = sprintf('%s in %s:%i', ...
            ME.message, ME.stack(1).file, ME.stack(1).line);
        else
          msg = ME.message;
        end
      else msg = res.Details.Output; end
      fprintf(1, '%10s %7s %s %s\n', ...
        res.Name, res.Details.Status, char(res.Details.Code), cleanupcomment(msg));
    end
  end
  disp('Totals:')
  fprintf(1, '  %i Passed, %i Failed, %i Incomplete.\n', ...
    sum([ r.Passed ]), sum([ r.Failed ]), sum([ r.Incomplete ]));
  fprintf(1, '  %f  seconds testing time.\n', sum([ r.Duration ]));

% ------------------------------------------------------------------------------
function [passed, failed, output] = runtest_sandbox(ln)
  passed=false; failed=false; output='';
  fprintf(1, '.');
  try
    clear ans
    output = evalc(ln);
    result = ans;
    if (isnumeric(result) || islogical(result)) && all(result(:))
      passed=true;
    elseif ischar(result) && any(strcmpi(result,{'OK','passed'}))
      passed=true;
    else
      failed=true;
    end
  catch ME
    output = ME;
    failed = true;
  end
    
  
