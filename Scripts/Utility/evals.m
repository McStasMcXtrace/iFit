function [passed, failed, output] = evals(ln)
% EVALS Evaluate expression similarly as EVALC but in a sandbox.
%   passed = EVALS(expr) evaluates expression, and returns 'true' when the result
%   is true, or generates 'ans' as non zero, or the output contains 'OK' or
%   'passed' as first word.
%
%   [passed, failed, output] = EVALS(expr) returns 'failed' as true when the
%   result was null/not OK, or failed execution. The output contains the
%   evaluation result, or a MException object (when failed execution).

  passed=false; failed=false; output='';
  if iscellstr(ln)
    ln = sprintf('%s; ', ln{:});
  end
  if ~ischar(ln)
    error([ mfilename ': expect a char/cellstr as expression to evaluate' ]);
  end

  try
    clear ans
    output = evalc(ln);
    if exist('ans','var')
      result = ans;
    else
      result = 1;
    end
    if (isnumeric(result) || islogical(result)) && all(result(:))
      passed=true;
    elseif ischar(result) && any(strcmpi(strtok(result),{'OK','passed'}))
      passed=true;
    else
      failed=true;
    end
  catch ME
    output = ME;
    failed = true;
  end
