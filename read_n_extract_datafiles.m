function[PSDValues]=read_n_extract_datafiles(Thresholds,filename1,filename2,filename3,varargin)
% This function is being called from scorematic.m
global EMG_SAMPLES EEG_SAMPLES EEG_TIMESTAMPS EPOCHSIZE State INDEX EEG_CHANNEL...
    EPOCHtime SampFreq1 FileFlag Statetime Statenum EPOCHnum Indices Fs EMG_CHANNEL stateTrack stateHistory
global D_lo D_hi T_lo T_hi S_lo S_hi B_lo B_hi
global EEG_Fc EMG_Fc
global EEG_highpass_enable EEG_HP_Fc EEG_Notch_enable EMG_Notch_enable EMG_LP_Fc EMG_lowpass_enable
global Input3_enable Input3_LP_enable Input3_LP_Fc Input3_Notch_enable Input3_HP_enable Input3_HP_Fc INPUT3_TIMESTAMPS INPUT3_SAMPLES INPUT3_CHANNEL
global Input4_enable Input4_HP_enable Input4_HP_Fc Input4_Notch_enable Input4_LP_enable Input4_LP_Fc INPUT4_TIMESTAMPS INPUT4_SAMPLES INPUT4_CHANNEL
global sampfactor exactLow exactHi boundIndex
% Take the variables out of the VARARGIN depending on its length.
waithandle=[];
waithandle= waitbar(0,'Loading the file ..... ');pause(0.5)

if length(varargin) == 3
    filename4=char(varargin(1,1));
    lowertimestamp = str2num(char(varargin(1,end-1)));
    uppertimestamp = str2num(char(varargin(1,end)));
%     lowerbound=str2num(char(varargin(1,end-1)));
%     upperbound=str2num(char(varargin(1,end)));
else
    lowertimestamp = str2num(char(varargin(1,end-1)));
    uppertimestamp = str2num(char(varargin(1,end)));
%     lowerbound=str2num(char(varargin(1,end-1)));
%     upperbound=str2num(char(varargin(1,end)));
    filename4=[];
end
handles=guihandles(scorematic);
cwd = pwd;
cd(tempdir);
%pack; %NOT COMMENTED IN ORIGINAL
cd(cwd);
% chkbox_value=get(handles.checkbox1,'Value');
nsamp = []; % Number of samples per bin of time
%ORIGINAL code at bottom of m-file.
if isequal(FileFlag, 2) % AD System
    try
        if EMG_CHANNEL <= 5
            [Timestamps,Samples,SF,Nsamp]=Crx2Mat(filename1,lowertimestamp+EMG_CHANNEL-1,...
                uppertimestamp+EMG_CHANNEL-1);
        else % Channel 6 have timestamps skipped by 1 from channel 5, hence the offset changes
            [Timestamps,Samples,SF,Nsamp]=Crx2Mat(filename1,lowertimestamp+EMG_CHANNEL,...
                uppertimestamp+EMG_CHANNEL);
        end
    catch %#ok<CTCH>
        fprintf('There is an error in Crx2Mat function. Check values passed into the function \n' );
        fprintf('\n Maybe the channel numbers do not match the actual channel number\n\n');
        rethrow(lasterror); %#ok<LERR>
    end
    waitbar(0.1,waithandle,'Converting EMG from CRextract format to Matlab dataformat ...');
    figure(waithandle),pause(0.2),
    samples=double(Samples(:)');
    clear Samples
    SampFreq1=SF;
    if boundIndex == 1
        DS = (1:1:10);
        DSampSF = SampFreq1./DS;
        indSampfactor = find(DSampSF >= 250);
        Fs = round(DSampSF(indSampfactor(end)));
        sampfactor = DS(indSampfactor(end));
        msgbox({['Orginal Sampling Rate:  ' num2str(SampFreq1) 'Hz'];...
            ['Down-Sampled Sampling Rate:  ' num2str(Fs) 'Hz']; ['Sampling Factor:  ' num2str(sampfactor) '']});
    end

    waitbar(0.4,waithandle,'Filtering the EMG data ...'); 
    figure(waithandle),pause(0.2),
    % Incase you are looking at stimulation file, then you dont have to
    % filter out the high frequency data. WE just use unfiltered data
    [path,name,ext]=fileparts(filename1);
    switch ext
        case '.stim'
            [EMG_TIMESTAMPS,EMG_SAMPLES]=Crxread(double(Timestamps),samples,sampfactor,Nsamp, exactLow, exactHi);
        case '.emg'
            [EMG_TIMESTAMPS,EMG_SAMPLES] = Crxread(double(Timestamps),samples,sampfactor,Nsamp, exactLow, exactHi);


            %  High pass filter for EMG signals
            [Bhigh,Ahigh]=ellip(7,1,60, EMG_Fc/(SampFreq1/2),'high');  % Default setting implements high pass filter with 30hz cutoff
            filtered_samples=(filter(Bhigh,Ahigh,EMG_SAMPLES))/4;
            %  OPTIONAL lowpass filter for EMG signals
            if EMG_lowpass_enable>0
                [EMG_Blow,EMG_Alow] = ellip(7,1,60, EMG_LP_Fc/(SampFreq1/2));   % Default is OFF
                filtered_samples = filter(EMG_Blow,EMG_Alow, filtered_samples);
            end
            % Optional EMG 60Hz Notch Filter
            if EMG_Notch_enable > 0
                woB = 60/(SampFreq1/2);
                [B_EMG_Notch,A_EMG_Notch] =  iirnotch(woB, woB/35);   % Default is OFF
                filtered_samples = filter(B_EMG_Notch,A_EMG_Notch, filtered_samples);
            end
            EMG_TIMESTAMPS = EMG_TIMESTAMPS(1:sampfactor:end);
            EMG_SAMPLES = filtered_samples(1:sampfactor:end);
    end
    
elseif isequal(FileFlag, 1) % Neuralynx System
    waitbar(0.1,waithandle,'Converting EMG from Neuralynx CSC to Matlab data ...');
    figure(waithandle),pause(0.2),
    [Timestamps,SF,Samples] = Nlx2MatCSC(filename1,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
    samples = double(Samples(:)');
    clear Samples
    SampFreq1 = SF(1);
    if boundIndex == 1
        DS = (1:1:10);
        DSampSF = SampFreq1./DS;
        indSampfactor = find(DSampSF >= 250);
        Fs = round(DSampSF(indSampfactor(end)));
        sampfactor = DS(indSampfactor(end));
         msgbox({['Orginal Sampling Rate:  ' num2str(SampFreq1) 'Hz'];...
            ['Down-Sampled Sampling Rate:  ' num2str(Fs) 'Hz']; ['Sampling Factor:  ' num2str(sampfactor) '']});
    end
    % Precise time stamps should be calculated here:
    [EMG_TIMESTAMPS,EMG_SAMPLES] = generate_timestamps_from_Ncsfiles(Timestamps, samples, exactLow, exactHi, nsamp);
    %Set up EMG filters:    
    %  High pass filter for EMG signals
    [Bhigh,Ahigh]=ellip(7,1,60, EMG_Fc/(SampFreq1/2),'high');  % Default setting implements high pass filter with 30hz cutoff
    waitbar(0.4,waithandle,'Filtering the EMG data ...'); 
    figure(waithandle),pause(0.2),
    physInput = 1;  %Needed to select proper error box in HeaderADBit.
    ADBit2uV = HeaderADBit(filename1, physInput);    %Calls a function to extract the AD Bit Value.
    EMG_SAMPLES = EMG_SAMPLES * ADBit2uV;   %Convert EMG amplitude of signal from AD Bits to microvolts.
    filtered_samples = filter(Bhigh,Ahigh,EMG_SAMPLES);
    % OPTIONAL lowpass filter for EMG signals
    if EMG_lowpass_enable>0
        [EMG_Blow,EMG_Alow] = ellip(7,1,60, EMG_LP_Fc/(SampFreq1/2));   % Default is OFF
        filtered_samples = filter(EMG_Blow,EMG_Alow, filtered_samples);
    end
    % Optional 60Hz Notch Filter
    if EMG_Notch_enable > 0
        woB = 60/(SampFreq1/2);
        [B_EMG_Notch,A_EMG_Notch] =  iirnotch(woB, woB/35);   % Default is OFF
        filtered_samples = filter(B_EMG_Notch,A_EMG_Notch, filtered_samples);
    end
    EMG_TIMESTAMPS = EMG_TIMESTAMPS(1:sampfactor:end);
    EMG_SAMPLES = filtered_samples(1:sampfactor:end);
    clear physInput ADBit2uV
elseif isequal(FileFlag, 3) % ASCII - Polysmith
    waitbar(0.1,waithandle,'Converting EMG from ASCII to Matlab data ...');
    figure(waithandle),pause(0.2),
    %Changed below from 200Hz
    SampFreq1 = 200;
    exactLow = (exactLow + 7200 * SampFreq1 * (boundIndex - 1));
    exactHi = (exactLow + exactHi - 1);
    SampA = dlmread(filename1);
    %Below is in Karin's version
    filenameB = strrep(filename1, '-1', '-2');
    SampB = dlmread(filenameB);
    Samples = SampA - SampB;
    clear SampA SampB
    Samples = Samples(exactLow:exactHi);
    samples = double(Samples(:));
    clear Samples
    physInput = 1;  %Needed to select proper error box in HeaderADBit.
    samples = AsciiPolysmithADBit(filename1, physInput, samples);
    
    tStep = 1/SampFreq1;
    tLow = exactLow/SampFreq1;
    tHi = exactHi/SampFreq1;
    Timestamps = (tLow:tStep:tHi);
    sampfactor = 1;
   
    Fs = round(SampFreq1/sampfactor);
    % Confirming that the downsampling is not too low
    if boundIndex == 1
        DS = (1:1:10);
        DSampSF = SampFreq1./DS;
        indSampfactor = find(DSampSF >= 200);
        Fs = round(DSampSF(indSampfactor(end)));
        sampfactor = DS(indSampfactor(end));
         msgbox({['Orginal Sampling Rate:  ' num2str(SampFreq1) 'Hz'];...
            ['Down-Sampled Sampling Rate:  ' num2str(Fs) 'Hz']; ['Sampling Factor:  ' num2str(sampfactor) '']});
    end

    % Precise time stamps should be calculated here:
    EMG_TIMESTAMPS = Timestamps;
    EMG_SAMPLES = samples;
    %Set up EMG filters    
    [Bhigh,Ahigh]=ellip(7,1,60, EMG_Fc/(SampFreq1/2),'high');  % Default setting implements high pass filter with 30hz cutoff
    %[Bhigh,Ahigh,Blow,Alow] = filterDefinition;
    waitbar(0.4,waithandle,'Filtering the EMG data ...'); 
    figure(waithandle),pause(0.2),
    filtered_samples = filter(Bhigh,Ahigh,EMG_SAMPLES);
    % OPTIONAL lowpass filter for EMG signals
    if EMG_lowpass_enable>0
        [EMG_Blow,EMG_Alow] = ellip(7,1,60, EMG_LP_Fc/(SampFreq1/2));   % Default is OFF
        filtered_samples = filter(EMG_Blow,EMG_Alow, filtered_samples);
    end
    % Optional 60Hz Notch Filter
    if EMG_Notch_enable > 0
        woB = 60/(SampFreq1/2);
        [B_EMG_Notch,A_EMG_Notch] =  iirnotch(woB, woB/35);   % Default is OFF
        filtered_samples = filter(B_EMG_Notch,A_EMG_Notch, filtered_samples);
    end
    EMG_TIMESTAMPS = EMG_TIMESTAMPS(1:sampfactor:end);
    EMG_SAMPLES = filtered_samples(1:sampfactor:end);
elseif isequal(FileFlag, 4) % ASCII - EMZA    
    waitbar(0.6,waithandle,' Converting EMG from ASCII format to Matlab dataformat ...');
    figure(waithandle),pause(0.2)
    SampFreq1 = 200;
    exactLow = (exactLow + 7200 * SampFreq1 * (boundIndex - 1));
    exactHi = (exactLow + exactHi - 1);

    % Code for extracting the data from the ASCII files:
    fileA = fopen(filename1);
    tempSampAtext = textscan(fileA, '%s',1);
    tempSampAnum = textscan(fileA, '%f %f', 'delimiter', '"');
    SampA = tempSampAnum{1,2}(:);
    clear tempSampAtext tempSampAnum
    fclose(fileA);

    filenameB = strrep(filename1, '-1', '-2');  % Automatically determines the name of the 2nd file of the pair.
    fileB = fopen(filenameB);
    tempSampBtext = textscan(fileB, '%s',1);
    tempSampBnum = textscan(fileB, '%f %f', 'delimiter', '"');
    SampB = tempSampBnum{1,2}(:);
    clear tempSampBtext tempSampBnum
    fclose(fileB);

    Samples = SampA - SampB;
    clear SampA SampB
    Samples = Samples(exactLow:exactHi);
    samples = double(Samples(:));
    clear Samples
    tStep = 1/SampFreq1;
    tLow = exactLow/SampFreq1;
    tHi = exactHi/SampFreq1;
    Timestamps = (tLow:tStep:tHi);
    sampfactor = 1;
   
    Fs = round(SampFreq1/sampfactor);
    % Confirming that the downsampling is not too low
    if boundIndex == 1
        DS = (1:1:10);
        DSampSF = SampFreq1./DS;
        indSampfactor = find(DSampSF >= 200);
        Fs = round(DSampSF(indSampfactor(end)));
        sampfactor = DS(indSampfactor(end));
         msgbox({['Orginal Sampling Rate:  ' num2str(SampFreq1) 'Hz'];...
            ['Down-Sampled Sampling Rate:  ' num2str(Fs) 'Hz']; ['Sampling Factor:  ' num2str(sampfactor) '']});
    end

    % Precise time stamps should be calculated here:
    EMG_TIMESTAMPS = Timestamps;
    EMG_SAMPLES = samples;
    %Set up EMG filters    
    [Bhigh,Ahigh]=ellip(7,1,60, EMG_Fc/(SampFreq1/2),'high');  % Default setting implements high pass filter with 30hz cutoff
    %[Bhigh,Ahigh,Blow,Alow] = filterDefinition;
    waitbar(0.4,waithandle,'Filtering the EMG data ...'); 
    figure(waithandle),pause(0.2),
    filtered_samples = filter(Bhigh,Ahigh,EMG_SAMPLES);
    % OPTIONAL lowpass filter for EMG signals
    if EMG_lowpass_enable>0
        [EMG_Blow,EMG_Alow] = ellip(7,1,60, EMG_LP_Fc/(SampFreq1/2));   % Default is OFF
        filtered_samples = filter(EMG_Blow,EMG_Alow, filtered_samples);
    end
    % Optional 60Hz Notch Filter
    if EMG_Notch_enable > 0
        woB = 60/(SampFreq1/2);
        [B_EMG_Notch,A_EMG_Notch] =  iirnotch(woB, woB/35);   % Default is OFF
        filtered_samples = filter(B_EMG_Notch,A_EMG_Notch, filtered_samples);
    end
    EMG_TIMESTAMPS = EMG_TIMESTAMPS(1:sampfactor:end);
    EMG_SAMPLES = filtered_samples(1:sampfactor:end);
elseif isequal(FileFlag, 5) % Plexon System
    channel = [];
    while isempty(channel)
        prompt={'Enter channel to be used for EMG:'};
        dlgTitle='EMG Channel Select';
        lineNo=1;
        answer = inputdlg(prompt,dlgTitle,lineNo);
        channel = str2double(answer{1,1});
        clear answer prompt dlgTitle lineNo
    end
    waitbar(0.1,waithandle,'Importing EMG data from Plexon PLX file ...');
    figure(waithandle),pause(0.2),
    [adfreq, ~, Timestamps, nsamp, Samples] = plx_ad_v(filename, channel);
    %[Timestamps,SF,Samples] = Nlx2MatCSC(filename,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
    samples = double(Samples(:)');
    clear Samples
    SampFreq1 = adfreq;
    if BoundIndex == 1
        DS = (1:1:10);
        DSampSF = SampFreq1./DS;
        indSampfactor = find(DSampSF >= 600);
        Fs = round(DSampSF(indSampfactor(end)));
        sampfactor = DS(indSampfactor(end));
         msgbox({['Orginal Sampling Rate:  ' num2str(SampFreq1) 'Hz'];...
            ['Down-Sampled Sampling Rate:  ' num2str(Fs) 'Hz']; ['Sampling Factor:  ' num2str(sampfactor) '']});
    end
    % Precise time stamps should be calculated here:
    [EMG_TIMESTAMPS,EMG_SAMPLES] = generate_timestamps_from_Ncsfiles(Timestamps, samples, exactLow, exactHi, nsamp);
    % 'generate_timestamps...' converts the time for Neuralynx files to
    % seconds.  The next equation is used to correct the time stamps for
    % PLX files.
    EMG_TIMESTAMPS = EMG_TIMESTAMPS * 1000000;
    %Set up EMG filters    
    %  High pass filter for EMG signals
    [Bhigh,Ahigh]=ellip(7,1,60, EMG_Fc/(SampFreq1/2),'high');  % Default setting implements high pass filter with 30hz cutoff
    waitbar(0.4,waithandle,'Filtering the EMG data ...'); 
    figure(waithandle),pause(0.2),
    filtered_samples = filter(Bhigh,Ahigh,EMG_SAMPLES);
    % OPTIONAL lowpass filter for EMG signals
    if EMG_lowpass_enable>0
        [EMG_Blow,EMG_Alow] = ellip(7,1,60, EMG_LP_Fc/(SampFreq1/2));   % Default is OFF
        filtered_samples = filter(EMG_Blow,EMG_Alow, filtered_samples);
    end
    % Optional 60Hz Notch Filter
    if EMG_Notch_enable > 0
        woB = 60/(SampFreq1/2);
        [B_EMG_Notch,A_EMG_Notch] =  iirnotch(woB, woB/35);   % Default is OFF
        filtered_samples = filter(B_EMG_Notch,A_EMG_Notch, filtered_samples);
    end
    EMG_TIMESTAMPS = EMG_TIMESTAMPS(1:sampfactor:end);
    EMG_SAMPLES = filtered_samples(1:sampfactor:end);    
end
clear filtered_samples samples adfreq

%  ******    EEG FILE extraction   *********
% filename=get(handles.eegfile,'TooltipString');
% lowertimestamp=LBounds(boundIndex);
% uppertimestamp=UBounds(boundIndex);
if isequal(FileFlag, 2) % AD System
    try
        if EEG_CHANNEL <= 5
            [Timestamps,Samples,SF,Nsamp]=Crx2Mat(filename2,lowertimestamp+EEG_CHANNEL-1,...
                uppertimestamp+EEG_CHANNEL-1);
        else % Channel 6 have timestamps skipped by 1 from channel 5, hence the offset changes
            [Timestamps,Samples,SF,Nsamp]=Crx2Mat(filename2,lowertimestamp+EEG_CHANNEL,...
                uppertimestamp+EEG_CHANNEL);
        end
    catch %#ok<CTCH>
        fprintf( 'There is an error in Crx2Mat function. Check values passed into the function \n' );
        rethrow(lasterror); %#ok<LERR>
    end
    waitbar(0.6,waithandle,' Converting EEG from Crextract format to Matlab dataformat ...');
    figure(waithandle),pause(0.2)
    samples=double(Samples(:)');
    clear Samples
    SampFreq1=SF;
    waitbar(0.8,waithandle,'Filtering the EEG data ...'); 
    figure(waithandle),pause(0.2),
    [EEG_TIMESTAMPS,EEG_SAMPLES]=Crxread(double(Timestamps),samples,sampfactor,Nsamp, exactLow, exactHi);
    [Blow,Alow]=ellip(7,1,60, EEG_Fc/(SampFreq2/2));           % Default setting implements low pass filter with 30hz cutoff
    filtered_samples=filter(Blow,Alow,EEG_SAMPLES/4);
    %  OPTIONAL highpass filter for EEG signals
    if EEG_highpass_enable>0
        [EEG_Bhi,EEG_Ahi] = ellip(7,1,60, EEG_HP_Fc/(SampFreq1/2),'high');   % Default is OFF
        filtered_samples = filter(EEG_Bhi,EEG_Ahi, filtered_samples);
    end
    %  OPTIONAL 60Hz Notch filter for EEG signals
    if EEG_Notch_enable > 0
        wo = 60/(SampFreq1/2);
        [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
        filtered_samples = filter(B_EEG_Notch,A_EEG_Notch, filtered_samples);
    end
    EEG_TIMESTAMPS = EEG_TIMESTAMPS(1:sampfactor:end);
    EEG_SAMPLES = filtered_samples(1:sampfactor:end);
elseif isequal(FileFlag, 1) % Neuralynx System
    waitbar(0.6,waithandle,' Converting EEG from Neuralynx CSC to Matlab data ...');
    figure(waithandle),pause(0.2)
    [Timestamps,SF,Samples]=Nlx2MatCSC(filename2,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
    samples=double(Samples(:)');
    clear Samples
    SampFreq2=SF(1);
    Fs=round(SampFreq2/sampfactor);
    waitbar(0.8,waithandle,'Filtering the EEG data ...'); 
    figure(waithandle),pause(0.2),
    [EEG_TIMESTAMPS,EEG_SAMPLES] = generate_timestamps_from_Ncsfiles(Timestamps,samples, exactLow, exactHi, nsamp);
    physInput = 2;  %Needed to select proper error box in HeaderADBit.
    ADBit2uV = HeaderADBit(filename2, physInput);    %Calls a function to extract the AD Bit Value.
    EEG_SAMPLES = EEG_SAMPLES * ADBit2uV;   %Convert EEG amplitude of signal from AD Bits to microvolts.
    %  Low pass filter for EEG signals
    [Blow,Alow]=ellip(7,1,60, EEG_Fc/(SampFreq2/2));           % Default setting implements low pass filter with 30hz cutoff
    filtered_samples=filter(Blow,Alow,EEG_SAMPLES);
    %  OPTIONAL highpass filter for EEG signals
    if EEG_highpass_enable>0
        [EEG_Bhi,EEG_Ahi] = ellip(7,1,60, EEG_HP_Fc/(SampFreq2/2),'high');   % Default is OFF
        filtered_samples = filter(EEG_Bhi,EEG_Ahi, filtered_samples);
    end
    %  OPTIONAL 60Hz Notch filter for EEG signals
    if EEG_Notch_enable > 0
        wo = 60/(SampFreq2/2);
        [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
        filtered_samples = filter(B_EEG_Notch,A_EEG_Notch, filtered_samples);
    end
    EEG_TIMESTAMPS = EEG_TIMESTAMPS(1:sampfactor:end);
    EEG_SAMPLES = filtered_samples(1:sampfactor:end);
    clear physInput ADBit2uV
elseif isequal(FileFlag, 3) % ASCII - Polysmith
    waitbar(0.6,waithandle,' Converting EEG from ASCII format to Matlab dataformat ...');
    figure(waithandle),pause(0.2)

    SampFreq1 = 200;
    SampA = dlmread(filename2);
    filenameB = strrep(filename2, '-1', '-2');
    SampB = dlmread(filenameB);
    Samples = SampA - SampB;
    clear SampA SampB
    Samples = Samples(exactLow:exactHi);
    samples = double(Samples(:));
    clear Samples
   
    physInput = 2;  %Needed to select proper error box in AsciiPolysmithADBit.
    samples = AsciiPolysmithADBit(filename2, physInput, samples);
    
    Fs = round(SampFreq1/sampfactor);
    waitbar(0.8,waithandle,'Filtering the EEG data ...'); 
    figure(waithandle),pause(0.2),
    EEG_TIMESTAMPS = Timestamps;
    EEG_SAMPLES = samples;
    [Blow,Alow]=ellip(7,1,60, EEG_Fc/(SampFreq1/2));
    filtered_samples=filter(Blow,Alow,EEG_SAMPLES);
    %  OPTIONAL highpass filter for EEG signals
    if EEG_highpass_enable>0
        [EEG_Bhi,EEG_Ahi] = ellip(7,1,60, EEG_HP_Fc/(SampFreq1/2),'high');   % Default is OFF
        filtered_samples = filter(EEG_Bhi,EEG_Ahi, filtered_samples);
    end
    %  OPTIONAL 60Hz Notch filter for EEG signals
    if EEG_Notch_enable > 0
        wo = 60/(SampFreq1/2);
        [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
        filtered_samples = filter(B_EEG_Notch,A_EEG_Notch, filtered_samples);
    end
    EEG_TIMESTAMPS = EEG_TIMESTAMPS(1:sampfactor:end);
    EEG_SAMPLES = filtered_samples(1:sampfactor:end);
    
elseif isequal(FileFlag, 4) % ASCII - EMZA    
    waitbar(0.6,waithandle,' Converting EEG from ASCII format to Matlab dataformat ...');
    figure(waithandle),pause(0.2)
    SampFreq1 = 200;
    
    % Code for extracting the data from the ASCII files:
    fileA = fopen(filename2);
    tempSampAtext = textscan(fileA, '%s',1);
    tempSampAnum = textscan(fileA, '%f %f', 'delimiter', '"');
    SampA = tempSampAnum{1,2}(:);
    clear tempSampAtext tempSampAnum
    fclose(fileA);

    filenameB = strrep(filename2, '-1', '-2');  % Automatically determines the name of the 2nd file of the pair.
    fileB = fopen(filenameB);
    tempSampBtext = textscan(fileB, '%s',1);
    tempSampBnum = textscan(fileB, '%f %f', 'delimiter', '"');
    SampB = tempSampBnum{1,2}(:);
    clear tempSampBtext tempSampBnum
    fclose(fileB);
    
    Samples = SampA - SampB;
    clear SampA SampB
    Samples = Samples(exactLow:exactHi);
    samples = double(Samples(:));
    clear Samples
   
    Fs = round(SampFreq1/sampfactor);
    waitbar(0.8,waithandle,'Filtering the EEG data ...'); 
    figure(waithandle),pause(0.2),
    EEG_TIMESTAMPS = Timestamps;
    EEG_SAMPLES = samples;
    [Blow,Alow]=ellip(7,1,60, EEG_Fc/(SampFreq1/2));
    filtered_samples=filter(Blow,Alow,EEG_SAMPLES);
    %  OPTIONAL highpass filter for EEG signals
    if EEG_highpass_enable>0
        [EEG_Bhi,EEG_Ahi] = ellip(7,1,60, EEG_HP_Fc/(SampFreq1/2),'high');   % Default is OFF
        filtered_samples = filter(EEG_Bhi,EEG_Ahi, filtered_samples);
    end
    %  OPTIONAL 60Hz Notch filter for EEG signals
    if EEG_Notch_enable > 0
        wo = 60/(SampFreq1/2);
        [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
        filtered_samples = filter(B_EEG_Notch,A_EEG_Notch, filtered_samples);
    end
    EEG_TIMESTAMPS = EEG_TIMESTAMPS(1:sampfactor:end);
    EEG_SAMPLES = filtered_samples(1:sampfactor:end);
elseif isequal(FileFlag, 5) % Plexon
    channel = [];
    while isempty(channel)
        prompt={'Enter channel to be used for EEG:'};
        dlgTitle='EEG Channel Select';
        lineNo=1;
        answer = inputdlg(prompt,dlgTitle,lineNo);
        channel = str2double(answer{1,1});
        clear answer prompt dlgTitle lineNo
    end
    waitbar(0.1,waithandle,'Importing EEG data from Plexon PLX file ...');
    figure(waithandle),pause(0.2),
    [adfreq, ~, Timestamps, nsamp, Samples] = plx_ad_v(filename, channel);
    %[Timestamps,SF,Samples] = Nlx2MatCSC(filename,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
    samples = double(Samples(:)');
    clear Samples
    SampFreq2 = adfreq;
    Fs=round(SampFreq2/sampfactor);
    waitbar(0.8,waithandle,'Filtering the EEG data ...'); 
    figure(waithandle),pause(0.2),
    [EEG_TIMESTAMPS,EEG_SAMPLES] = generate_timestamps_from_Ncsfiles(Timestamps,samples, exactLow, exactHi, nsamp);
    % 'generate_timestamps...' converts the time for Neuralynx files to
    % seconds.  The next equation is used to correct the time stamps for
    % PLX files.
    EEG_TIMESTAMPS = EEG_TIMESTAMPS * 1000000;
    %  Low pass filter for EEG signals
    [Blow,Alow]=ellip(7,1,60, EEG_Fc/(SampFreq2/2));           % Default setting implements low pass filter with 30hz cutoff
    filtered_samples=filter(Blow,Alow,EEG_SAMPLES);
    %  OPTIONAL highpass filter for EEG signals
    if EEG_highpass_enable>0
        [EEG_Bhi,EEG_Ahi] = ellip(7,1,60, EEG_HP_Fc/(SampFreq2/2),'high');   % Default is OFF
        filtered_samples = filter(EEG_Bhi,EEG_Ahi, filtered_samples);
    end
    %  OPTIONAL 60Hz Notch filter for EEG signals
    if EEG_Notch_enable > 0
        wo = 60/(SampFreq2/2);
        [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
        filtered_samples = filter(B_EEG_Notch,A_EEG_Notch, filtered_samples);
    end
    EEG_TIMESTAMPS = EEG_TIMESTAMPS(1:sampfactor:end);
    EEG_SAMPLES = filtered_samples(1:sampfactor:end);
end
clear filtered_samples samples SF
%freqVersusTimePowerMap(EEG_SAMPLES, Fs, EEG_Fc);
%*******************************************************************

EPOCHtime=[]; EPOCHnum=[];
if isempty(filename4)==0    %Checks to see if Training file exists
    try
        [t_stamps]=xlsread(filename4);  %Imports Training file data
    catch
        uiwait(errordlg('Check if the file is saved in Microsoft Excel format.',...
            'ERROR','modal'));
    end
    EPOCHtime = t_stamps(1:end,2);  % Makes column of downsampled time steps new vector 
    EPOCHnum = t_stamps(1:end,3);   % Makes column of scored states for each time step above
    clear t_stamps
end
fprintf('          Calculating the FFT and getting the power values....\n');
% This is to find the EPOCHSIZE of 10 sec
index=find((EMG_TIMESTAMPS(1)+9.999 < EMG_TIMESTAMPS) & (EMG_TIMESTAMPS < EMG_TIMESTAMPS(1)+10.001));
if(isempty(index)) == 1
    index=find((EMG_TIMESTAMPS(1)+9.99 < EMG_TIMESTAMPS) & (EMG_TIMESTAMPS < EMG_TIMESTAMPS(1)+10.01));
end
diff= EMG_TIMESTAMPS(index(1):index(end)) - (EMG_TIMESTAMPS(1)+10);
[minimum,ind]=min(abs(diff));
try
    EPOCHSIZE=index(ind);
catch
        fprintf('There is an error in calculating the EPOCHSIZE of 10sec in read_n_extract_datafiles\n');
end
%Take the fft of the entire timeframe we have, calcuate power in bins
[PSDValues]=fft_psd_and_statescore_of_epoch(Thresholds);
waitbar(1,waithandle, 'Finished converting.. Now Loading the data ..');
figure(waithandle), pause(0.3);
close(waithandle);