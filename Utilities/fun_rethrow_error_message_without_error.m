function fun_rethrow_error_message_without_error(ME, hyperlinkOnQ)
if nargin < 2
    hyperlinkOnQ = true;
end
if hyperlinkOnQ
    disp(getReport(ME, 'extended', 'hyperlinks', 'on'));
else
    disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
end
end