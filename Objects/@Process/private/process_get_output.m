function out = process_get_output(stream)
  % process_get_output: get the stream (InputStream) content[private]
  
    out = '';
    available = 0;
    try
      available = stream.available;
    end
    
    % return when nothing to read or invalid stream
    if available <= 0, return; end
    
    % Read the content of the stream.
    %
    % EPIC FAIL 1:
    % The following method would be nice, but fails:
    %   target_buffer = javaArray('java.lang.Byte', available);
    %   stream.read(target_buffer, 0, available);
    %
    % Indeed, as matlab converts any Java array into a Matlab class, 
    % the read method can not be used with additional parameters (signature not found)
    % see http://www.mathworks.com/matlabcentral/answers/66227-syntax-for-call-to-java-library-function-with-byte-reference-parameter
    %
    % EPIC FAIL 2:
    % using readLine from a BufferedReader stalls after a few iterations
    % Reader       = java.io.InputStreamReader(stream);
    % bufferReader = java.io.BufferedReader(Reader);
    % readLine(bufferReader)
    
    % we use a for loop to read all bytes, one by one (sadly).
    out = zeros(1, available);
    for index=1:available
      out(index) = read(stream);
    end
    
    out = strrep(char(out), sprintf('\n\r'), sprintf('\n'));

  end
