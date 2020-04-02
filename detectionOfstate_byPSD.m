% #####################################################################
function[nChangedIndices] = detectionOfstate_byPSD(operateOn,nIndexofepochs,...
    detectState,Thresholds)
% This function is being called from fft_psd_and_statescore_of_epoch.m file
%The function accepts Indices of the epochs which are potentially a
%particular epochs, and calculates 1sec FFT for the 10sec long epochs and 
%finds out what the majority of the epoch power distribution looks like. 
%Depending on that it assigns the state. It returns the indices of epochs 
%which would be then named as the changed states.

global Fs EEG_SAMPLES EMG_SAMPLES EPOCHSIZE Statetime State EEG_TIMESTAMPS averageSTpower StdDevforAverageSTpower

nIterations=length(nIndexofepochs);
nChangedIndices=[];
Start_pt=[];
if operateOn == 'eeg'    
    loopindex=1;
    while loopindex <= nIterations
        originalIndex=nIndexofepochs(loopindex);
        EPOCH_startpoint = (EPOCHSIZE * (originalIndex-1)) + 1;
        
        if (double(EEG_TIMESTAMPS(EPOCH_startpoint)) + 9.9) < double(EEG_TIMESTAMPS(end))
            maxind=9;
        else
            Epochsize=double(EEG_TIMESTAMPS(end)) - double(EEG_TIMESTAMPS(EPOCH_startpoint));
            maxind=floor(Epochsize/10);
        end
        ind=0; 
        for ind=0:maxind
            unfiltered_data=double(EEG_SAMPLES(EPOCH_startpoint+round(ind*EPOCHSIZE/10)...
                : EPOCH_startpoint+round((ind+1)*EPOCHSIZE/10)));
            W=ones(size(unfiltered_data)); W=W(1:end-2);
            %W=hann(length(unfiltered_data)-2);
            df=length(unfiltered_data);
            [Pxx2,F2]=spectrum(double(unfiltered_data),df,0,W(1:end-1),Fs);
            indexdelta=[];indextheta=[];indexsigma=[];
            %For the EEG signal
            indexdelta=find(F2(1)+0.4<F2 & F2 < F2(1)+4);   % delta band 0.4 -4.9 Hz
            indextheta=find(F2(1)+5 < F2 & F2 < F2(1)+9);  % theta band 5-10 Hz
            indexsigma=find(F2(1)+10< F2 & F2 < F2(1)+14);  % sigma band 11-14 Hz
            deltapower(ind+1)=sum(Pxx2(indexdelta))/df *2;
            thetapower(ind+1)=sum(Pxx2(indextheta))/df *2;
            sigmapower(ind+1)=sum(Pxx2(indexsigma))/df *2;
        end
        switch detectState
            case 'tr'
                sigmaTimesTheta = sigmapower .* thetapower;
                nEpochsOfHighSigma=find(sigmaTimesTheta > (5*averageSTpower));% + 10*StdDevforAverageSTpower));
                if length(nEpochsOfHighSigma) >= 6
                    nChangedIndices = [nChangedIndices; originalIndex];
                end
                %Original published logic:
%                 nEpochsOfHighSigma=find(sigmapower > Thresholds.sigma2SD);
%                 if length(nEpochsOfHighSigma) >= 5
%                     nChangedIndices = [nChangedIndices; originalIndex];
%                 end
        end
        loopindex=loopindex+1;
    end
else        % operateOn =='emg'
    loopindex=1;
    switch detectState
        case 'qw'
            emgthreshold = 2.5*Thresholds.emg;
            % 3.5*EMG threshold was based on heureistic measure
        case 'aw'
            emgthreshold = Thresholds.emg;
    end
    while loopindex <= nIterations
        originalIndex=nIndexofepochs(loopindex);
        EPOCH_startpoint = (EPOCHSIZE * (originalIndex-1)) + 1;
%         if (loopindex ==432) 
%             EPOCH_startpoint
%             EEG_TIMESTAMPS(EPOCH_startpoint)
%             EEG_TIMESTAMPS(end)
%             if (double(EEG_TIMESTAMPS(EPOCH_startpoint)) + 9.9) < double(EEG_TIMESTAMPS(end))
%                 maxind=9
%             else
%                 Epochsize=double(EEG_TIMESTAMPS(end)) - double(EEG_TIMESTAMPS(EPOCH_startpoint))
%                 maxind=floor(Epochsize/10)
%             end
%             
%             emg_power=zeros(1,9);
%             whos EEG_SAMPLES EMG_SAMPLES           
%             for ind=0:maxind
%                 EPOCH_startpoint+round((ind+1)*EPOCHSIZE/10)
%                 emg_data=double(EMG_SAMPLES(EPOCH_startpoint+round(ind*EPOCHSIZE/10)...
%                     : EPOCH_startpoint+round((ind+1)*EPOCHSIZE/10)));
%                 absolute_emgdata=abs(emg_data).^2;             % absolute square
%                 emg_power(ind+1)=sum(absolute_emgdata)/length(absolute_emgdata); % sum of all squared Vj's
%                 ind
%             end
%             
%             nEpochsOfHighEmg=find(emg_power > emgthreshold);
%             if length(nEpochsOfHighEmg) >= ceil(maxind/2); % majority of 1 sec epochs within 10sec
%                 nChangedIndices=[nChangedIndices;originalIndex]
%             end
%             loopindex=loopindex+1;
%             
%         end
        %Start_pt=[Start_pt; EEG_TIMESTAMPS(EPOCH_startpoint)];
        ind=0; 
        if EPOCH_startpoint <= length(EEG_TIMESTAMPS)
            if (double(EEG_TIMESTAMPS(EPOCH_startpoint)) + 9.9) < double(EEG_TIMESTAMPS(end));
                maxind=9;
            else
                Epochsize=double(EEG_TIMESTAMPS(end)) - double(EEG_TIMESTAMPS(EPOCH_startpoint));
                maxind=floor(Epochsize/10);
            end
        end
            if (maxind < 1)
                loopindex = loopindex +1;
            else
                emg_power=zeros(1,9);
                for ind=0:maxind
                    st_sec = EPOCH_startpoint+round(ind*EPOCHSIZE/10);
                    end_sec =EPOCH_startpoint+round((ind+1)*EPOCHSIZE/10);
                    if (end_sec > length(EMG_SAMPLES))
                        end_sec = length(EMG_SAMPLES);
                    end
                    emg_data=double(EMG_SAMPLES(st_sec:end_sec));
    %                 emg_data=double(EMG_SAMPLES(EPOCH_startpoint+round(ind*EPOCHSIZE/10)...
    %                     : EPOCH_startpoint+round((ind+1)*EPOCHSIZE/10)));
                    absolute_emgdata=abs(emg_data).^2;             % absolute square
                    emg_power(ind+1)=sum(absolute_emgdata)/length(absolute_emgdata); % sum of all squared Vj's
                end

                nEpochsOfHighEmg=find(emg_power > emgthreshold);
                if length(nEpochsOfHighEmg) >= ceil(maxind/2); % majority of 1 sec epochs within 10sec
                    nChangedIndices=[nChangedIndices;originalIndex];
                end
                loopindex=loopindex+1;
            end
        
    end
end     
fprintf(' Finished with detection\n');