% This function draws a line from the first temporal point to the last in
% the variable XYF. z is a vector that contains the missing points that
% will be filled in with the line of best fit.
% XYF - current set of points
% z - frame locations for missing points
% New_XYF - XYF plus points @ frame=z
function [New_XYF] = CTv0_Func_FillMissingPoints(XYF,z)
    Avg = mean(XYF);
    m = XYF(end,:) - XYF(1,:); % slope calculated from first to last point of Temp
    m  = m/m(3);
    % ------------------------------------------------
    T = (z-Avg(3))/m(3); % t parameter at frame locations
    Missing_XYF = Avg + T'*m; % Calc missing coords of line of best fit
    
    XYF = [XYF; Missing_XYF];
    [~,idx]=sort(XYF(:,3)); 
    New_XYF=round(XYF(idx,:)); 
end

