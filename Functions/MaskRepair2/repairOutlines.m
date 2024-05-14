function repairOutlines(D,firstLast,opts)

% initiate loop variables
inter = 0; % variable to interupt the loop upon hitting case 's' = stop
N = 500; % number of points for the final outline; change if needed.
changed = 1;

ff = firstLast(1);
while ff <= firstLast(2)
% for ff=1:1:length(D.outlines) % flip through the outlines 1 by 1.
    % Import image and outline data
    I = importFrame(D.images,ff); % Import frame of interest corresponding to the outline
    outline = importdata(quickfilepath(D.outlines(ff))); % Import the outline of interest
    outlinePrev = outline; % store a copy of the outline as a backup in case the changes need to be deleted.

    % Initiate while loop variables
    inloop = 1; % While in the loop

    while inloop
                close all
                % Plot the mask and outline
                plotFrameWithOutline(I,outline,opts)
                title([sprintf('Outline %i',ff),"Press 'a' to repair, 'p' to see previous, 'z' to undo, 'return' to save & move to next, 's' to stop"],'interpreter','none')

            keyChar = waitForKeyPress();
            switch keyChar
                case 'a' % Replace points
                    % outlinePrev = outline;
                    outline = fixOutline(I,outline,opts);
                    % ff = ff+1;
                    changed = 1;

                case 'z' % Undo
                    outline = outlinePrev; % Get initial, unchanged outline
                    changed = 0;

                case 's' % Stop
                    inter = 1; % Interrupt outer loop
                    break;

                case 'p' % use outline from previous frame
                    ff = ff-1;
                    changed = 0;
                    break
                                    
                otherwise 
                    % ff = ff+1;
                    changed = 1;
                    break; % Leave inner loop and proceed to saving the outline
            
            end
     
    end

    if inter
        close all
        disp('Stopping')
        break
    end

    if changed == 1
        if size(outline,1)>1 % if the outline file is empty because no mask was created, then simply keep moving forward
            % Save outlines with evenly spaced points
            outline = curvspace(outline,N);
        else% if the outline file is empty because no mask was created, then simply resave the file with 0,0 in the first line and keep moving forward
            outline = [0,0];
        end
    % export individual outline to a .csv file (also re-export those that did
    % not need to be fixed and that were skipped when simply hitting return.
    filenameOutline = sprintf('iface_%05g', ff); % Number sequentially
    writematrix(outline,fullfile(D.outlines(1).folder,[filenameOutline,'.csv']))
        
    ff = ff+1;
    end 
end

% Clear the workspace when all done
clf
close all
fprintf('All masks repaired\n');

end