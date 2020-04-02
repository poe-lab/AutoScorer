% ###################################################################
function[PowerSpectralDensityValues]=fft_psd_and_statescore_of_epoch(Thresholds)
% This function does the FFT and does all the power calculations here.
% This function is being called from the M file - read_n_extract_datafiles.m
% this was updated last 3-17-04 because of a problem with the autoscorer
% lines 15, 279-281, 301-303 added and emg stuff added to lines 14, 313, 314, and 316
% on 4-24-06 by Theresa in order to get mean EMG values in excel file of the PSDvalues
% output 

% This version has detection of QW done before RE and TR detection
global EMG_SAMPLES EEG_SAMPLES EEG_TIMESTAMPS EPOCHSIZE State INDEX Statetime...
    EPOCHtime SampFreq1 Fs Statenum EPOCHnum Indices Flag stateTrack stateHistory
persistent SquaresumSigma SumSigma LengthSigma SquaresumTheta SumTheta...
    LengthTheta SquaresumDelta SumDelta LengthDelta SquaresumBeta SumBeta LengthBeta SquaresumEMG SumEMG LengthEMG
global D_lo D_hi T_lo T_hi S_lo S_hi B_lo B_hi averageSTpower StdDevforAverageSTpower

%NOT COMMENTED IN ORIGINAL
handles=guihandles(scorematic);
all_arraylength=ceil((double(EEG_TIMESTAMPS(end))-double(EEG_TIMESTAMPS(1)))/10);
Statenum=repmat(1,all_arraylength,1); Statetime=repmat(1,all_arraylength,1);
State=[repmat('A',all_arraylength,1) repmat('W',all_arraylength,1)];
if Flag==1
    stateTrack(1:all_arraylength) = {'B1'};
    stateHistory(1:all_arraylength) = {'B1'};
end
st_power=repmat(1,all_arraylength,1);dt_ratio=repmat(1,all_arraylength,1);
P_emg=repmat(1,all_arraylength,1);P_delta=repmat(1,all_arraylength,1);
P_theta=repmat(1,all_arraylength,1);P_sigma=repmat(1,all_arraylength,1);
P_beta=repmat(1,all_arraylength,1);
Indices=[];
EPOCH_StartPoint=1;INDEX=1;
nEpochPoints=length(EEG_TIMESTAMPS);
df=EPOCHSIZE;

while EPOCH_StartPoint <= nEpochPoints
    st_pt=EPOCH_StartPoint;
    end_pt=st_pt+EPOCHSIZE-1;  
    if length(EEG_TIMESTAMPS) - end_pt <= 0
        % If the current epoch is smaller then 10 sec, then dont score at
        % all... Just break out of the loop.. ( Changed on 06/20/03
        break;
    elseif length(EEG_TIMESTAMPS) - end_pt <= EPOCHSIZE/9
        end_pt = length(EEG_TIMESTAMPS);
    end
    Train_index=find(double(EEG_TIMESTAMPS(st_pt))-3 <= EPOCHtime & EPOCHtime <= double(EEG_TIMESTAMPS(st_pt))+3);
    if isempty(Train_index)==0
        Indices=[Indices;Train_index INDEX];
    end
    % This is taking integral in time domain  
    Vj=double(EMG_SAMPLES(st_pt:end_pt));
    absVj=abs(Vj).^2;                       % absolute square
    P_emg(INDEX)=sum(absVj)/length(absVj);  % sum of all squared Vj's
    
    % This is calculating power in frequency domain
    fft_in=double(EEG_SAMPLES(st_pt:end_pt));
    windowsize =length(fft_in);
    if df < windowsize
        df = windowsize;
    end
    
%     h = spectrum.welch('Hann', ones(windowsize,1), 0);  % Form: h = spectrum.welch('Hann',window,100*noverlap/window);
%     hpsd = psd(h, double(fft_in), 'NFFT', df, 'Fs', Fs);    %Form: hpsd = psd(h,x,'NFFT',nfft,'Fs',Fs);
%     Pxx2 = hpsd.Data;
%     F2 = hpsd.Frequencies;
    % SPECTRUM is obsolete. Replaced with above 4 lines that does the equivalent.
    [Pxx2,F2]=spectrum(double(fft_in),df,0,ones(windowsize,1),Fs);
    % ******  [P,F] = SPECTRUM(X,NFFT,NOVERLAP,WINDOW,Fs)
    
    %For the EEG signal
    index_delta=[];index_theta=[];index_sigma=[]; index_beta=[];
    index_delta=find(F2(1)+D_lo< F2 & F2 < F2(1)+D_hi);      % Default delta band 0.4 -4 Hz
    index_theta=find(F2(1)+T_lo< F2 & F2 < F2(1)+T_hi);    % Default theta band 5-9 Hz
    index_sigma=find(F2(1)+S_lo< F2 & F2 < F2(1)+S_hi);     % Default sigma band 10-14 Hz
    index_beta =find(F2(1)+B_lo< F2 & F2 < F2(1)+B_hi);     % Default Beta band 15-20 Hz
    
    P_delta(INDEX)=sum(Pxx2(index_delta))/df *2;  
    P_theta(INDEX)=sum(Pxx2(index_theta))/df *2;
    P_sigma(INDEX)=sum(Pxx2(index_sigma))/df *2;
    P_beta(INDEX)=sum(Pxx2(index_beta))/df *2;
    
    st_power(INDEX)=abs(P_sigma(INDEX)*P_theta(INDEX));   % Used to indicate waking
    dt_ratio(INDEX)=abs(P_delta(INDEX)/P_theta(INDEX));   
    
    Statetime(INDEX)=EEG_TIMESTAMPS(EPOCH_StartPoint);
    EPOCH_StartPoint=st_pt+EPOCHSIZE;        % Keep this at the end of this function
    INDEX=INDEX+1;
    warning('off','MATLAB:divideByZero');
end
% New Parameter on next line:
averageSTpower = mean(st_power);
StdDevforAverageSTpower = std(st_power);

if Flag==1 %The program is in auto-scoring mode.
  
    % Score epochs either as REM or QS depending on their DT ratio-EMG
    State(1:INDEX-1,1)='Q';State(1:INDEX-1,2)='S';
    Statenum(1:INDEX-1)=2;
    stateTrack(1:INDEX-1)={'Q1'};
    stateHistory(1:INDEX-1)={'Q1'};
    
    index_dt_below_thresh=find(dt_ratio < Thresholds.deltatheta & P_emg < Thresholds.emg);
    if isempty(index_dt_below_thresh) == 0
        State(index_dt_below_thresh,1)='R';State(index_dt_below_thresh,2)='E';
        Statenum(index_dt_below_thresh)=3;
        stateTrack(index_dt_below_thresh)={'R1'};
        for i=1:length(index_dt_below_thresh)
            stateHistory{index_dt_below_thresh(i)} = [stateHistory{index_dt_below_thresh(i)} 'R1'];
        end
    end
    
    % Now seperate out points which can be coined as WAKING
    index_all_waking=find(P_emg > Thresholds.emg);
    % Now from all those epochs which can be termed as QW/ AW based on
    % STthresh
    index_st=find(st_power(index_all_waking) < Thresholds.sigmatheta);
    index_st_below_thresh=index_all_waking(index_st);
    if isempty(index_st_below_thresh) == 0
        State(index_st_below_thresh,1)='A';State(index_st_below_thresh,2)='W';
        Statenum(index_st_below_thresh)=1;
        stateTrack(index_st_below_thresh)={'A1'};
        for i=1:length(index_st_below_thresh)
            stateHistory{index_st_below_thresh(i)} = [stateHistory{index_st_below_thresh(i)} 'A1'];
        end
    end
    
%     index_q=find(dt_ratio(index_st_below_thresh) < Thresholds.deltatheta);
%     index_quietwaking=index_st_below_thresh(index_q);
%     if isempty(index_quietwaking) == 0
%         State(index_quietwaking,1)='Q';State(index_quietwaking,2)='W';
%         Statenum(index_quietwaking)=4;
%     end
    
    % Check if the threshold is given & then look for the UNhooked epochs  
    if isempty(Thresholds.unhooked)==0
        index_unhooked=find(P_emg < 1 & st_power < Thresholds.unhooked);  % ### Change value accordingly
        if isempty(index_unhooked) == 0
            State(index_unhooked,1)='U';State(index_unhooked,2)='H';
            Statenum(index_unhooked)=5;
            stateTrack(index_unhooked)={'U1'};
            for i=1:length(index_unhooked)
                stateHistory{index_unhooked(i)}= [stateHistory{index_unhooked(i)} 'U1'];
            end
        end
    end
%     fprintf('calculated the first round states \n');
    
    % For absolute detection of AW states from states with very high EMG as
    % well as high Sigma * Theta value which are scored as QS before this
%     fprintf('**1**\n  ');
    index_qs=find(Statenum == 2);
    if isempty(index_qs) == 0
        high_st_emg=find(P_emg(index_qs) > Thresholds.emg & ...
            st_power(index_qs) > Thresholds.sigmatheta);
        if isempty(high_st_emg) == 0
            index_should_be_aw = index_qs(high_st_emg);
            
            % Now using these indices, we have to run FFT on 1sec epochs to
            % find out if majority of the 10sec epoch has AW or QS patterns
            operateOn='emg';  % On what should the 1sec FFT be done on...
            detectState='aw';
            [indices_of_aw] = detectionOfstate_byPSD(operateOn,index_should_be_aw,...
                detectState,Thresholds);
            qs_aw=length(indices_of_aw);
            State(indices_of_aw,1)='A';State(indices_of_aw,2)='W';
            Statenum(indices_of_aw)=1;
            stateTrack(indices_of_aw)={'A2'};
            for i=1:length(indices_of_aw)
                stateHistory{indices_of_aw(i)}= [stateHistory{indices_of_aw(i)} 'A2'];
            end
            %From those which are not AW, check if they have low D/T power
            % and hence to be called as REM 
            index_not_aw = ismember(index_should_be_aw,indices_of_aw);
            % This gives epochs which are not AW among those from which
            % could be AW
            index_could_be_rem= find(index_not_aw == 0);
            index_lowdt = find(dt_ratio(index_could_be_rem) < Thresholds.deltatheta);
            index_is_rem = index_could_be_rem(index_lowdt);
            ind_rem=find(P_emg(index_is_rem) < (Thresholds.emg + 0.1*Thresholds.emg));
            index_is_rem=index_is_rem(ind_rem);
            State(index_is_rem,1)='R';State(index_is_rem,2)='E';
            Statenum(index_is_rem)=3;
            stateTrack(index_is_rem)={'R2'};
            for i=1:length(index_is_rem)
                stateHistory{index_is_rem(i)}= [stateHistory{index_is_rem(i)} 'R2'];
            end
        end
    end
            
%     fprintf('**3** \n  ');   
    % To see if there are any Transition to REM state within QS states
    index_qs=find(Statenum == 2);
    if isempty(index_qs) == 0
        
% This is the modified logic that was being used until 3/6/2013. Original published logic has been restored.        
%         if index_qs(1)== 1
%             index_qs=index_qs(2:end);
%         end
%         
%        %Run thru the - detectionOfstate_byPSD for TR detection
%         operateOn='eeg';
%         detectState='tr';
%         [indices_of_tr]=detectionOfstate_byPSD(operateOn,...
%             index_qs,detectState,Thresholds);
%         qs_tr = length(indices_of_tr);
%         State(indices_of_tr,1)='T';State(indices_of_tr,2)='R';
%         Statenum(indices_of_tr)=6;
%         stateTrack(indices_of_tr)={'T1'};
%         for i=1:length(indices_of_tr)
%             stateHistory{indices_of_tr(i)}= [stateHistory{indices_of_tr(i)} 'T1'];
%         end
                
            
        
      %ORIGINAL Published logic:  
        if index_qs(1)== 1
            index_qs=index_qs(2:end);
        end
        ind_sig_gr=find(P_sigma(index_qs) >Thresholds.sigma3SD);
        if isempty(ind_sig_gr) == 0
            index_poss_tr = index_qs(ind_sig_gr);
            ind_prev_qs_re=find(Statenum((index_poss_tr)-1) == 2 | Statenum((index_poss_tr)-1) == 3);  % Previous state is QS or RE
            if isempty(ind_prev_qs_re) == 0
                index_transition_rem=index_poss_tr(ind_prev_qs_re);
                if isempty(index_transition_rem)== 0
                   %Run thru the - detectionOfstate_byPSD for TR detection
                    operateOn='eeg';
                    detectState='tr';
                    [indices_of_tr]=detectionOfstate_byPSD(operateOn,...
                        index_transition_rem,detectState,Thresholds);
                    qs_tr=length(indices_of_tr);
                    State(indices_of_tr,1)='T';State(indices_of_tr,2)='R';
                    Statenum(indices_of_tr)=6;
                    stateTrack(indices_of_tr)={'T1'};
                    for i=1:length(indices_of_tr)
                        stateHistory{indices_of_tr(i)}= [stateHistory{indices_of_tr(i)} 'T1'];
                    end
                end
            end
        end
    end
    
    %fprintf('**4** \n  ');
    % To see if there is any QW within the AW or QS states by checking
    % their sigma, delta and theta levels. They should be below the
    % threshold set by the user
    index_quietwaking =[];
    indexaw=find(Statenum == 1);   indexqs=find(Statenum == 2);
    if isempty(indexqs)==0
        smallDeltaQSIndx = find(P_delta(indexqs) < Thresholds.delta);
        indexqs = indexqs(smallDeltaQSIndx);
    end
    index_aw_qs=[indexaw ;indexqs];
    index_aw_qs=sort(index_aw_qs);
    
    if isempty(index_aw_qs)==0
        low_sig_del_theta=find(P_theta(index_aw_qs) < Thresholds.theta & P_sigma(index_aw_qs) < Thresholds.sigma);
        if isempty(low_sig_del_theta) == 0
            index_quietwaking=index_aw_qs(low_sig_del_theta);
            aw_qs_qw=length(index_quietwaking);
            %Run thru the - detectionOfstate_byPSD for AW detection
            operateOn='emg';
            detectState='qw';
            [indices_of_aw]=detectionOfstate_byPSD(operateOn,...
                index_quietwaking,detectState,Thresholds);
            index_not_aw=setxor(index_quietwaking,indices_of_aw);
            fprintf(' Passed point 4\n');
            % Gives indices which are not AW
            State(index_not_aw,1)='Q';State(index_not_aw,2)='W';
            Statenum(index_not_aw)=4;
            stateTrack(index_not_aw)={'W1'};
            for i=1:length(index_not_aw)
                stateHistory{index_not_aw(i)}= [stateHistory{index_not_aw(i)} 'W1'];
            end
        end
    end
    
    fprintf('**2**\n  ');
    % For absolute REM detection with the REM states
    index_rem=find(Statenum == 3);        
    if isempty(index_rem)==0
        if (index_rem(1)==1 || index_rem(1)==2)
            index_rem=index_rem(2:end);
            if isempty(index_rem)==0
                if index_rem(1)== 2
                    index_rem=index_rem(2:end);
                end
            end
        end
        if isempty(index_rem)==0
            if index_rem(1)== 3
                index_rem=index_rem(2:end);
            end
        end
        
        if isempty(index_rem)==0
            A= find(Statenum(index_rem-1)==1 | Statenum(index_rem-1)==4);
            B= find(Statenum(index_rem-2)==1 | Statenum(index_rem-2)==4);
            C= find(Statenum(index_rem-3)==1 | Statenum(index_rem-3)==4);
            D = intersect(A, B);
            ind_prev_3waking = intersect(C, D);
            index_is_qw = index_rem(ind_prev_3waking);
            State(index_is_qw,1)='Q';State(index_is_qw,2)='W';
            Statenum(index_is_qw)=4;
            stateTrack(index_is_qw)={'W2'};
            for i=1:length(index_is_qw)
                stateHistory{index_is_qw(i)}= [stateHistory{index_is_qw(i)} 'W2'];
            end
        end
    end

    index_QW = find(Statenum == 4);
    if isempty(index_QW)==0
        hi_theta=find(P_theta(index_QW) > Thresholds.theta);
         if isempty(hi_theta) == 0
            index_intermedWaking = index_QW(hi_theta);
            State(index_intermedWaking,1) = 'I'; State(index_intermedWaking,2) = 'W';
            Statenum(index_intermedWaking) = 8;
            stateTrack(index_intermedWaking)={'I1'};
            for i=1:length(index_intermedWaking)
                stateHistory{index_intermedWaking(i)}= [stateHistory{index_intermedWaking(i)} 'I1'];
            end
         end
    end
        
end   % End for the if Flag==1 

 fprintf(' Passed point 5\n');
% Since the size of these arrays are pre-sized, we just check if array
% length and the INDEX number matches... if not, then resize it correctly

if INDEX-1 < all_arraylength
    
    Statenum=Statenum(1:INDEX-1); Statetime=Statetime(1:INDEX-1);
    State=State(1:INDEX-1,:);  
    st_power=st_power(1:INDEX-1);dt_ratio=dt_ratio(1:INDEX-1);
    P_emg=P_emg(1:INDEX-1);P_delta=P_delta(1:INDEX-1);
    P_theta=P_theta(1:INDEX-1);P_sigma=P_sigma(1:INDEX-1);
    P_beta=P_beta(1:INDEX-1);
    %P_alpha=P_alpha(1:INDEX-1);
end
if isempty(SquaresumSigma) == 0
    
    SquaresumSigma=SquaresumSigma + sum(P_sigma.^2);
    SumSigma=SumSigma + sum(P_sigma);
    LengthSigma=LengthSigma + length(P_sigma);
    
    SquaresumDelta=SquaresumDelta + sum(P_delta.^2);
    SumDelta=SumDelta + sum(P_delta);
    LengthDelta=LengthDelta + length(P_delta);
    
    SquaresumTheta=SquaresumTheta + sum(P_theta.^2);
    SumTheta=SumTheta + sum(P_theta);
    LengthTheta=LengthTheta + length(P_theta);
    
    SquaresumBeta=SquaresumBeta + sum(P_beta.^2);
    SumBeta=SumBeta + sum(P_beta);
    LengthBeta=LengthBeta + length(P_beta);
    
    SquaresumEMG=SquaresumEMG + sum(P_emg.^2);
    SumEMG=SumEMG + sum(P_emg);
    LengthEMG=LengthEMG + length(P_emg);
         
    
else
    SquaresumSigma = sum(P_sigma.^2);
    SumSigma = sum(P_sigma);
    LengthSigma = length(P_sigma);
    
    SquaresumDelta = sum(P_delta.^2);
    SumDelta = sum(P_delta);
    LengthDelta = length(P_delta);
    
    SquaresumTheta = sum(P_theta.^2);
    SumTheta = sum(P_theta);
    LengthTheta = length(P_theta);
    
    SquaresumBeta = sum(P_beta.^2);
    SumBeta = sum(P_beta);
    LengthBeta = length(P_beta);
    
    SquaresumEMG = sum(P_emg.^2);
    SumEMG = sum(P_emg);
    LengthEMG = length(P_emg);
    
    
    
end

% Create a structure to be passed out of this function having all the power
% levels of different bands of interest

SquareSum=struct('delta',SquaresumDelta, 'theta',SquaresumTheta,...
    'sigma',SquaresumSigma, 'beta',SquaresumBeta, 'emg',SquaresumEMG );
Sum=struct('delta',SumDelta, 'theta',SumTheta, 'sigma',SumSigma, 'beta',SumBeta, 'emg',SumEMG);
TotalPoints=struct('delta',LengthDelta, 'theta',LengthTheta,...
    'sigma',LengthSigma, 'beta',LengthBeta, 'emg',LengthEMG);

PowerSpectralDensityValues = struct( 'emg',P_emg, 'sigmatheta',st_power,...
    'deltatheta',dt_ratio, 'delta',P_delta, 'theta',P_theta,'sigma',P_sigma, ...
     'beta',P_beta, 'squaresum',SquareSum,'sum',Sum, 'nPoints',TotalPoints);
    % 'beta',P_beta,'alpha',P_alpha, 'squaresum',SquareSum,'sum',Sum, 'nPoints',TotalPoints);
fprintf('Finished fft_psd file\n');