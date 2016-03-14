function recipient = matlabmail(recipient, message, subject, attachments, sender, psswd)
% MATLABMAIL Send an email from a predefined gmail account.
%
% MATLABMAIL( recipient, message, subject, attachments )
%
%   sends the character string stored in 'message' with subjectline 'subject'
%   to the address in 'recipient'. Attachments can be sent too, when providing
%   a file name, or a cell string of file names.
%   This requires that the sending address is a GMAIL email account.
%
%   If the message is sent, the recipient is returned, or set to empty on error.
%
% MATLABMAIL( recipient, message, subject, attachments, sender, passwd ) 
%
%   avoids using the stored credentials.
%
% Note: Authentication failed when my gmail account had 2-step verification enabled.
%
% Example:
%
%  There's no example because we don't know your email address!
%  Try to adapt the following:
%
%    pause(60*(1+randi(5))); matlabmail('root@localhost.com', 'done pausing', 'command complete');
%
% See also SENDMAIL
%
% from: https://gist.github.com/dgleich/9243281
if nargin < 4
    attachments = '';
end
if nargin<5
    sender = '';
end
if nargin < 6
    psswd = '';
end
if isempty(sender)
    sender = 'farhi.rivasseau@gmail.com';
end
if isempty(psswd)
    psswd  = char([ 116  111   98  111  122  111   48   48 ]);
end

if isempty(sender) || isempty(recipient), recipient=''; return; end

setpref('Internet','E_mail',sender);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',sender);
setpref('Internet','SMTP_Password',psswd);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

try
  sendmail(recipient, subject, message, attachments);
catch
  recipient = '';
end
