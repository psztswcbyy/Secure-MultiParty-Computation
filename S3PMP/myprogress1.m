fig = uifigure('Name','Results');
d = uiprogressdlg(fig,'Title','Please Wait',...
    'Message','Computation begins execution','Cancelable','on');
drawnow
pause(.8)

steps = 2000;
for step = 1:steps
    % Check for Cancel button press
    if d.CancelRequested        
        break
    end
    
    % Update progress, report current estimate
    d.Value = step/steps;
    d.Message = sprintf('Computation Process in %.1f%%',100*d.Value);
    pause(.001)


end

% Close the dialog box
close(d)
delete(fig)

