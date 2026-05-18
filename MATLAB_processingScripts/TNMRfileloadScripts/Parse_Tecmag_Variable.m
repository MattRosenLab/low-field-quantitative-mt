function val=Parse_Tecmag_Variable(lineval)
% Parse_Tecmag_Variable.m: Takes a character string lineval, which may
% include a letter at the end indicating units, and converts to a value of 
% type double, returned as newval. If lineval is already of type double, it
% passes through untouched.
% -- Dave Korenchan, 12 Apr 2026

% Extract from cell, if input is a cell
if iscell(lineval)
    lineval=lineval{:};
end

% Remove any whitespace/control characters
lineval=lineval(find(~isstrprop(lineval,'cntrl')&~isstrprop(lineval,'wspace')));

% ID units and convert to appropriate double value
if ischar(lineval)
    if contains(lineval,'u')
        val=str2double(strtok(lineval,'u'))*1e-6;
    elseif contains(lineval,'m')
        val=str2double(strtok(lineval,'m'))*1e-3;
    elseif contains(lineval,'s')
        val=str2double(strtok(lineval,'s'));
    else
        val=str2double(lineval);
    end  
else
    val=lineval;
end