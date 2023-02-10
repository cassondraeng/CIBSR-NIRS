function [rawSplit] = extract_task_period(raw)
% This function extracts the during task period only 


rawSplit = [];
for i =1:numel(raw)
    try 
        if (ismember('qr_during', raw(i).stimulus.keys))
            tStart = raw(i).stimulus.qr_during.onset;
            tEnd = raw(i).stimulus.qr_during.onset + raw(i).stimulus.qr_during.dur;
            idx = find((raw(i).time >= tStart) & raw(i).time < tEnd); % find index of time vector between start and end
            time = raw(i).time(idx);  
            data = raw(i).data(idx,:);
            
            rawQR = raw(i);
            rawQR.time = time;
            rawQR.data = data;
            rawQR.demographics.task = 'qr_during';

            rawSplit = [rawSplit rawQR]; % add to dataset 

        elseif (ismember('in_during', raw(i).stimulus.keys))
            tStart = raw(i).stimulus.in_during.onset;
            tEnd = raw(i).stimulus.in_during.onset + raw(i).stimulus.in_during.dur;
            idx = find((raw(i).time >= tStart) & raw(i).time < tEnd); % find index of time vector between start and end
            time = raw(i).time(idx);  
            data = raw(i).data(idx,:);
            
            rawIN = raw(i);
            rawIN.time = time;
            rawIN.data = data;
            rawIN.demographics.task = 'in_during';

            rawSplit = [rawSplit rawIN]; % add to dataset 
        elseif (ismember('out_during', raw(i).stimulus.keys))
            tStart = raw(i).stimulus.out_during.onset;
            tEnd = raw(i).stimulus.out_during.onset + raw(i).stimulus.out_during.dur;
            idx = find((raw(i).time >= tStart) & raw(i).time < tEnd); % find index of time vector between start and end
            time = raw(i).time(idx);  
            data = raw(i).data(idx,:);
            
            rawOUT = raw(i);
            rawOUT.time = time;
            rawOUT.data = data;
            rawOUT.demographics.task = 'out_during';

            rawSplit = [rawSplit rawOUT]; % add to dataset 
        else
            fprintf('%s, no during in file\n', raw(i).demographics.subject);
        end

    catch e
        disp(raw(i).demographics.subject); disp('caused error')
        disp(i);
       % list_error{kk} = i;
       % kk = kk+1;
        disp(e.identifier);
        disp(e.message);
        disp('=========================================================')
        continue
    end
end

end