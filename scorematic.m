%##########################################################################
%                              Auto-Scorer 2008
%
%                                 Created by 
%                               Brooks A. Gross
%
%                            SLEEP AND MEMORY LAB
%                           UNIVERSITY OF MICHIGAN
%##########################################################################
% DESCRIPTION:
% This GUI is used to automatically score sleep states for EEG and EMG 
% data.  It takes in data from the Neuralynx or AD system files in 2 hour 
% blocks.  The filter settings file is automatically loaded in the 
% pre-scored mode when a training file is present and updates the filter 
% settings GUI components and variables.  The auto-scored states are 
% automatically saved into an excel file.
%##########################################################################
% VERSION 2
% Modified by Brooks A. Gross on Jan-23-2008
% --Filter information now in a separate file that is automatically named 
%   based on UI Training file name. It is automatically imported to
%   scorematic.m during pre-scored mode.
% --Added all new filter options in Sleep Scorer to the scorematic GUI.
% --Added pop-up menu for setting downsampling value.
% --Added 'Set/Reset' button to load in all default values. This button
%   MUST be pressed before proceeding with any task in GUI.
% --'Read in AD files' checkbox now defaults as off since Neuralynx is used
%   more often.
%
% VERSION 1
% Modified by Theresa Bjorness on April-24-2006
% --Lines 344 and 404 added and 343 and 403 altered (\n changed to \t) to 
%   get mean EMG in the PSDvalues output excel spreadsheet.
%##########################################################################
function varargout = scorematic(varargin)
% SCOREMATIC Application M-file for scorematic.fig
%    FIG = SCOREMATIC launch scorematic GUI.
%    SCOREMATIC('callback_unhooked', ...) invoke the named callback.
% Last Modified by GUIDE v2.5 25-Mar-2010 09:58:59
% 
if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse',varargin{:});
    
    % Generate a structure of handles to pass to callbacks, and store it.  
    handles = guihandles(fig);
    guidata(fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end
    
end
% #########################################################################
%                          Code for Opening Files
% #########################################################################
function emgfile_open_Callback(hObject, eventdata, handles)
% This function is used to call the EMG file name.
working_dir=pwd;
current_dir='C:\SleepData\Datafiles';
cd(current_dir);
[filename, pathname] = uigetfile({'*.dat;*.Ncs;*.ncs; *.emg; *.ascii',...
        'All data files (*.dat, *.Ncs, *.ncs, *.emg, *.ascii)'},'Pick an EMG data file');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    emgfile= fullfile(pathname, filename);
    set(handles.emgfile,'string',filename);
    set(handles.emgfile,'Tooltipstring',emgfile);
end
% -------------------------------------------------------------------------
function eegfile_open_Callback(hObject, eventdata, handles)
% This function is used to call the EEG file name.
working_dir=pwd;
current_dir='C:\SleepData\Datafiles';
cd(current_dir);
[filename, pathname] = uigetfile({'*.dat;*.Ncs;*.ncs;*.eeg; *.ascii',...
        'All data files (*.dat, *.Ncs, *.ncs, *.eeg, *.ascii)'}, 'Pick an EEG data file');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    eegfile= fullfile(pathname, filename);
    set(handles.eegfile,'string',filename);
    set(handles.eegfile,'Tooltipstring',eegfile); 
end
% -------------------------------------------------------------------------
function trainingfile_open_Callback(hObject, eventdata, handles)
% This function is used to call the training file name.
working_dir=pwd;
current_dir='C:\SleepData\Results';
cd(current_dir);
[filename, pathname] = uigetfile('*.xls', 'Pick the training file you created');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    trainingfile= fullfile(pathname, filename);
    set(handles.trainingfile,'string',filename);
    set(handles.trainingfile,'Tooltipstring',trainingfile);
end
% -------------------------------------------------------------------------
function tstampsfile_open_Callback(hObject, eventdata, handles)
% This function is used to call the time stamp file name.
working_dir=pwd;
current_dir='C:\SleepData\Timestampfiles';
cd(current_dir);
[filename, pathname] = uigetfile('*.xls', 'Pick the timestamp file for these datafiles');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    timestampfile= fullfile(pathname, filename);
    set(handles.tstampsfile,'string',filename);
    set(handles.tstampsfile,'Tooltipstring',timestampfile);
end
% -------------------------------------------------------------------------
function autoscoredfile_open_Callback(hObject, eventdata, handles)
% This function is used to call the autoscored file name.
working_dir=pwd;
current_dir='C:\SleepData\Results';
cd(current_dir);
[filename, pathname] = uigetfile('*.xls', 'Pick a Autoscored file');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    autoscoredfile= fullfile(pathname, filename);
    set(handles.autoscoredfile,'string',filename);
    set(handles.autoscoredfile,'Tooltipstring',autoscoredfile);
end
% --------------------------------------------------------------------
function input3file_open_Callback(hObject, eventdata, handles)
working_dir=pwd;
current_dir='C:\SleepData\Datafiles';
cd(current_dir);
[filename, pathname] = uigetfile({'*.dat;*.Ncs;*.ncs;*.eeg; *.ascii',...
        'All data files (*.dat, *.Ncs, *.ncs, *.eeg, *.ascii)'}, 'Pick a file for Input 3');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    Input3file= fullfile(pathname, filename);
    set(handles.Input3file,'string',filename);
    set(handles.Input3file,'Tooltipstring',Input3file); 
end
% --------------------------------------------------------------------
function input4file_open_Callback(hObject, eventdata, handles)
working_dir=pwd;
current_dir='C:\SleepData\Datafiles';
cd(current_dir);
[filename, pathname] = uigetfile({'*.dat;*.Ncs;*.ncs;*.eeg; *.ascii',...
        'All data files (*.dat, *.Ncs, *.ncs, *.eeg, *.ascii)'}, 'Pick a file for Input 4');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    Input4file= fullfile(pathname, filename);
    set(handles.Input4file,'string',filename);
    set(handles.Input4file,'Tooltipstring',Input4file); 
end
% #########################################################################
% function filemenu_Callback(hObject, eventdata, handles)
% #########################################################################
function varargout = fileSelectionMenu_Callback(h, eventdata, handles, varargin)
global FileFlag
switch get(handles.fileSelectionMenu,'Value')   
    case 1
        FileFlag = 0;
        errordlg('You must select a file type.','Error');
    case 2 % Neuralynx
        FileFlag = 1;
    case 3 % AD System
        FileFlag = 2;
    case 4 % ASCII - Polysmith
        FileFlag = 3;
    case 5 % ASCII - EMZA
        FileFlag = 4;
end

% #########################################################################
function Scorematic_ResizeFcn(hObject, eventdata, handles)
% #########################################################################
% Code below is for selecting the EEG and EMG channels.
% #########################################################################
function eeg_cb1_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EEG_CHANNEL = 1;
SetEegChannelsOff
set(gcbo,'Checked','ON');
%
function eeg_cb2_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EEG_CHANNEL = 2;
SetEegChannelsOff
set(gcbo,'Checked','ON');
%
function eeg_cb3_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EEG_CHANNEL = 3;
SetEegChannelsOff
set(gcbo,'Checked','ON');
%
function eeg_cb4_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EEG_CHANNEL = 4;
SetEegChannelsOff
set(gcbo,'Checked','ON');
%
function eeg_cb5_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EEG_CHANNEL = 5;
SetEegChannelsOff
set(gcbo,'Checked','ON');
%
function emg_cb1_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EMG_CHANNEL = 1;
SetEmgChannelsOff
set(gcbo,'Checked','ON');
%
function emg_cb2_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EMG_CHANNEL = 2;
SetEmgChannelsOff
set(gcbo,'Checked','ON');
%
function emg_cb3_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EMG_CHANNEL = 3;
SetEmgChannelsOff
set(gcbo,'Checked','ON');
%
function emg_cb4_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EMG_CHANNEL = 4;
SetEmgChannelsOff
set(gcbo,'Checked','ON');
%
function emg_cb5_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EMG_CHANNEL = 5;
SetEmgChannelsOff
set(gcbo,'Checked','ON');
%
function eeg_cb6_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EEG_CHANNEL = 6;
SetEegChannelsOff
set(gcbo,'Checked','ON');
%
function eeg_cb7_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EEG_CHANNEL = 7;
SetEegChannelsOff
set(gcbo,'Checked','ON');
%
function eeg_cb8_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EEG_CHANNEL = 8;
SetEegChannelsOff
set(gcbo,'Checked','ON');
%
function emg_cb6_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EMG_CHANNEL = 6;
SetEmgChannelsOff
set(gcbo,'Checked','ON');
%
function emg_cb7_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EMG_CHANNEL = 7;
SetEmgChannelsOff
set(gcbo,'Checked','ON');
%
function emg_cb8_Callback(hObject, eventdata, handles)
global EMG_CHANNEL EEG_CHANNEL 
EMG_CHANNEL = 8;
SetEmgChannelsOff
set(gcbo,'Checked','ON');
%
function SetEegChannelsOff()
handles=guihandles(scorematic);
set(handles.eeg_cb1,'Checked','OFF');set(handles.eeg_cb2,'Checked','OFF');
set(handles.eeg_cb3,'Checked','OFF');set(handles.eeg_cb4,'Checked','OFF');
set(handles.eeg_cb5,'Checked','OFF');set(handles.eeg_cb6,'Checked','OFF');
set(handles.eeg_cb7,'Checked','OFF');set(handles.eeg_cb8,'Checked','OFF');
%
function SetEmgChannelsOff()
handles=guihandles(scorematic);
set(handles.emg_cb1,'Checked','OFF');set(handles.emg_cb2,'Checked','OFF');
set(handles.emg_cb3,'Checked','OFF');set(handles.emg_cb4,'Checked','OFF');
set(handles.emg_cb5,'Checked','OFF');set(handles.emg_cb6,'Checked','OFF');
set(handles.emg_cb7,'Checked','OFF');set(handles.emg_cb8,'Checked','OFF');
%##########################################################################
%                      End of channel selection code
%##########################################################################
function prescored_Callback(hObject, eventdata, handles)

global EMG_SAMPLES EEG_SAMPLES EEG_TIMESTAMPS EPOCHSIZE State INDEX EPOCHtime...
    Flag Statetime Statenum EPOCHnum Indices Cmap X1 Y1 Z1 X2 Y2 Z2
global D_lo D_hi T_lo T_hi S_lo S_hi B_lo B_hi
global EEG_Fc EMG_Fc
global EEG_highpass_enable EEG_HP_Fc EEG_Notch_enable EMG_Notch_enable EMG_LP_Fc EMG_lowpass_enable
global Input3_enable Input3_LP_enable Input3_LP_Fc Input3_Notch_enable Input3_HP_enable Input3_HP_Fc INPUT3_TIMESTAMPS INPUT3_SAMPLES INPUT3_CHANNEL
global Input4_enable Input4_HP_enable Input4_HP_Fc Input4_Notch_enable Input4_LP_enable Input4_LP_Fc INPUT4_TIMESTAMPS INPUT4_SAMPLES INPUT4_CHANNEL
global exactLow exactHi boundIndex
% This file uses PERSISTENT variables, and hence we clear them on every new use
clear('fft_psd_and_statescore_of_epoch')
% Get UI file names

filename1=get(handles.emgfile,'TooltipString');
filename2=get(handles.eegfile,'TooltipString');
filename3=get(handles.tstampsfile,'TooltipString');
filename4=get(handles.trainingfile,'TooltipString');
working_dir=pwd;
current_dir='C:\SleepData\Results';
cd(current_dir);

if isempty(filename4) == 0 % If training filename is provided
    filename5 = regexprep(filename4, '.xls', '');
    headrFiltrs = xlsread([filename5 'Header.xls'],'b1:y1')'
    cd(working_dir);
    if headrFiltrs(1) ~= -1
        EMG_Fc = headrFiltrs(1);
        set(handles.EMG_cutoff, 'String', EMG_Fc);
    end
    if headrFiltrs(3) ~= -1
        EMG_LP_Fc = headrFiltrs(3);
        set(handles.EMG_Lowpass_cutoff, 'String', EMG_LP_Fc);
        EMG_lowpass_enable = 1;
        set(handles.EMG_LP_checkbox, 'Value',1);
    end
    if headrFiltrs(5) ~= -1
        EMG_Notch_enable=1;
        set(handles.Notch60EMG, 'Value',1);
    end
    if headrFiltrs(7) ~= -1
        EEG_Fc = headrFiltrs(7);
        set(handles.EEG_cutoff, 'String', EEG_Fc);
    end
    if headrFiltrs(9) ~= -1
        EEG_HP_Fc = headrFiltrs(9);
        set(handles.EEG_Highpass_cutoff, 'String', EEG_HP_Fc);
        EEG_highpass_enable=1;
        set(handles.EEG_Highpass_checkbox, 'Value',1);
    end
    if headrFiltrs(11) ~= -1
        EEG_Notch_enable=1;
        set(handles.Notch60EEG, 'Value',1);
    end
    if (headrFiltrs(13)~= -1||headrFiltrs(15)~= -1||headrFiltrs(17)~= -1)
        Input3_enable = 1;
        set(handles.Input3_checkbox, 'Value',1);
    end
    if headrFiltrs(13) ~= -1
        Input3_HP_Fc = headrFiltrs(13);
        set(handles.Input3_Highpass_cutoff, 'String', Input3_HP_Fc);
        Input3_HP_enable = 1;
        set(handles.Input3_HP_checkbox, 'Value',1);
    end
    if headrFiltrs(15) ~= -1
        Input3_LP_Fc = headrFiltrs(15);
        set(handles.Input3_Lowpass_cutoff, 'String', Input3_LP_Fc);
        Input3_LP_enable = 1;
        set(handles.Input3_LP_checkbox, 'Value',1);
    end
    if headrFiltrs(17) ~= -1
        Input3_Notch_enable = 1;
        set(handles.Notch60Input3, 'Value',1);
    end
    if (headrFiltrs(19)~= -1||headrFiltrs(21)~= -1||headrFiltrs(23)~= -1)
        Input4_enable = 1;
        set(handles.Input4_checkbox, 'Value',1);
    end
    if headrFiltrs(19) ~= -1
        Input4_HP_Fc = headrFiltrs(19);
        set(handles.Input4_Highpass_cutoff, 'String', Input4_HP_Fc);
        Input4_HP_enable = 1;
        set(handles.Input4_HP_checkbox, 'Value',1);
    end
    if headrFiltrs(21) ~= -1
        Input4_LP_Fc = headrFiltrs(21);
        set(handles.Input4_Lowpass_cutoff, 'String', Input4_LP_Fc);
        Input4_LP_enable = 1;
        set(handles.Input4_LP_checkbox, 'Value',1);
    end
    if headrFiltrs(23) ~= -1
        Input4_Notch_enable = 1;
        set(handles.Notch60Input4, 'Value',1);
    end
else
    cd(working_dir);
end
ST_thresh =900;  
DT_thresh =1.1;  
EMG_thresh=15;
% Add to the search path of Matlab, the path where it can find the saved
% workspace variables. Remove the path at the end of the program.

figure(1)
% Setup the colormap for the states.
Cmap(1,:)=[1 0.8 0];  % Yellow   => Active Waking
Cmap(2,:)=[0 0 1];    % Blue     => Quiet Sleep
Cmap(3,:)=[1 0 0];    % Red      => REM
Cmap(4,:)=[0 1 0.1];  % Green    => Quiet Waking
Cmap(5,:)=[0 0 0];    % Black    => Unhooked
Cmap(6,:)=[0 1 1];    % Cyan     => Trans REM
Cmap(7,:)= [0.85 0.85 0.85];    % Grey    => Cleared State
Cmap(8,:)=[1 1 1];    % White    => Intermediate Waking

% Setup the colormap for the preview stage of a 3-D plot.
Pmap(1,:)=[0.73 0.45 0.8]; Pmap(2,:)=[0.97 0.5 0.7]; Pmap(3,:)=[0.622 0.796 0.012];
Pmap(4,:)=[0.88 0.4 0.2]; Pmap(5,:)=[0.7 0.3 0.9]; Pmap(6,:)=[0.2 0.5 1];  
Pmap(7,:)=[0 1 0.5]; Pmap(8,:)=[0.98 0.23 0.573];

Flag=0; % Indicates that fft_psd_and_statescore_of_epoch.m is called for plotting
        % If Flag=1, it means its in auto scoring mode.
working_dir=pwd;
current_dir='C:\SleepData\Results';
cd(current_dir);
[filename,pathname] = uiputfile('*PSDvalues.xls','Save PSDvalues file as:');
cd(working_dir);
[path,name1,ext]=fileparts(filename1); % Breaks up the file name string into parts
figname=strcat(pathname,'sdt_TR_',name1); %Concatenates the strings in () into a new string
newname=strcat(pathname,filename);

fid=fopen(newname,'w');
fprintf(fid,'The mean + 2*Std.dev value of Sigma power is:\t');
fprintf(fid,'\n');
fprintf(fid,'Timestamp\t');
fprintf(fid,'P_delta\t');
fprintf(fid,'P_theta\t');
fprintf(fid,'P_sigma\t');
% fprintf(fid,'P_alpha\t');
fprintf(fid,'P_beta\t');
fprintf(fid,'P_emg\n');
fclose(fid);

% Read in the timestampsfilebutton file.
[tbounds] = xlsread(filename3);
lbound = tbounds(1:end,1);  ubound = tbounds(1:end,2);
exactLowIndx = tbounds(1:end,3); exactHiIndx = tbounds(1:end,4);
figure(1);     
whitebg(figure(1),[1,1,1])     %[.165,.465,.665])
set(gcf,'color',[187/255, 201/255, 214/255],'position',[672 533 560 420])
set(gca,'Xcolor',[0, 0, 0],'YColor',[0 0 0],'ZColor',[0 0 0])

figure(2);     
whitebg(figure(2),[1,1,1])     %[.165,.465,.665])
set(gcf,'color',[187/255, 201/255, 214/255],'position',[672 33 560 420])
set(gca,'Xcolor',[0, 0, 0],'YColor',[0 0 0],'ZColor',[0 0 0])

Sum_P_sigma=0; Length_P_sigma=0; Squaresum_P_sigma=0;
Mean_sigma=[]; Std_dev_sigma=[];

for boundIndex=1:length(ubound)
    lowerbound=lbound(boundIndex);
    upperbound=ubound(boundIndex);
    exactLow = exactLowIndx(boundIndex);
    exactHi = exactHiIndx(boundIndex);
    % Thresholds is made '[]' because its only initialized when we
    % autoscore it.. Its a structure having all threshold values
    Thresholds=[];
    if isempty(filename4) ==0 % If training filename is provided
        [PSDvalues] = read_n_extract_datafiles(Thresholds,filename1,filename2,filename3,filename4,num2str(lowerbound),num2str(upperbound));
    else                    % carry on without the file
        [PSDvalues] = read_n_extract_datafiles(Thresholds,filename1,filename2,filename3,num2str(lowerbound),num2str(upperbound));
    end
    %clear('read_n_extract_datafiles')
    S=[];C=[];
    lowdt=[];highdt=[];lowst=[];highst=[];lowemg=[];highemg=[];dtless=[];emghigh=[];
    S=ones(INDEX-1,1);
    S(:)=10;
    C=ones(length(S),3);
    C = C * 0.85;
    
%     lowdt=find(PSDvalues.deltatheta<DT_thresh & PSDvalues.sigmatheta<ST_thresh & PSDvalues.emg<EMG_thresh);
%     C(lowdt,1)=Pmap(1,1);C(lowdt,2)=Pmap(1,2);C(lowdt,3)=Pmap(1,3);
%     highdt=find(PSDvalues.deltatheta>DT_thresh & PSDvalues.sigmatheta<ST_thresh & PSDvalues.emg<EMG_thresh);
%     C(highdt,1)=Pmap(2,1);C(highdt,2)=Pmap(2,2);C(highdt,3)=Pmap(2,3);
%     lowst=find(PSDvalues.deltatheta>DT_thresh & PSDvalues.sigmatheta<ST_thresh & PSDvalues.emg>EMG_thresh);
%     C(lowst,1)=Pmap(3,1);C(lowst,2)=Pmap(3,2);C(lowst,3)=Pmap(3,3);
%     highst=find(PSDvalues.deltatheta>DT_thresh & PSDvalues.sigmatheta>ST_thresh & PSDvalues.emg>EMG_thresh);
%     C(highst,1)=Pmap(4,1);C(highst,2)=Pmap(4,2);C(highst,3)=Pmap(4,3);
%     lowemg=find(PSDvalues.deltatheta<DT_thresh & PSDvalues.sigmatheta>ST_thresh & PSDvalues.emg<EMG_thresh);
%     C(lowemg,1)=Pmap(5,1);C(lowemg,2)=Pmap(5,2);C(lowemg,3)=Pmap(5,3);
%     highemg=find(PSDvalues.deltatheta>DT_thresh & PSDvalues.sigmatheta>ST_thresh & PSDvalues.emg<EMG_thresh);
%     C(highemg,1)=Pmap(6,1);C(highemg,2)=Pmap(6,2);C(highemg,3)=Pmap(6,3);
%     dtless=find(PSDvalues.deltatheta<DT_thresh & PSDvalues.sigmatheta<ST_thresh & PSDvalues.emg>EMG_thresh);
%     C(dtless,1)=Pmap(7,1);C(dtless,2)=Pmap(7,2);C(dtless,3)=Pmap(7,3);
%     emghigh=find(PSDvalues.deltatheta<DT_thresh & PSDvalues.sigmatheta>ST_thresh & PSDvalues.emg>EMG_thresh);
%     C(emghigh,1)=Pmap(8,1);C(emghigh,2)=Pmap(8,2);C(emghigh,3)=Pmap(8,3);

    if isempty(Indices)==0
        S(Indices(:,2))=50; I=Indices(:,2);
        num=EPOCHnum(Indices(:,1));
        C(I,1)=Cmap(num,1);C(I,2)=Cmap(num,2);C(I,3)=Cmap(num,3);
    end
    i=1;
    fid=fopen(newname,'a+');
    for i=1:length(State)
        fprintf(fid,'%9.3f\t',Statetime(i));
        fprintf(fid,'%6.3f\t',PSDvalues.delta(i));
        fprintf(fid,'%6.3f\t',PSDvalues.theta(i));
        fprintf(fid,'%6.3f\t',PSDvalues.sigma(i));
        % fprintf(fid,'%6.3f\t',PSDvalues.alpha(i));
        fprintf(fid,'%6.3f\t',PSDvalues.beta(i));
        fprintf(fid,'%6.3f\n',PSDvalues.emg(i));
    end
    fclose(fid);
    figure(1),scatter3(PSDvalues.deltatheta,PSDvalues.sigmatheta,PSDvalues.emg,S,C,'o','filled');
    xlabel('Delta/Theta Power');ylabel('Sigma*Theta Power');zlabel('EMG Power');
    set(get(gca,'Title'),'Color','k')
    set(get(gca,'XLabel'),'Color','k')
    set(get(gca,'YLabel'),'Color','k')
    set(get(gca,'ZLabel'),'Color','k')
    set(gca,'Xcolor','k','YColor','k','ZColor','k')
    view(-65,22); hold on, pause(2)
    figure(2),scatter3(PSDvalues.sigma,PSDvalues.theta,PSDvalues.delta,S,C,'filled');
    xlabel('Sigma Power');ylabel('Theta Power');zlabel('Delta Power');
    set(get(gca,'Title'),'Color','k')
    set(get(gca,'XLabel'),'Color','k')
    set(get(gca,'YLabel'),'Color','k')
    set(get(gca,'ZLabel'),'Color','k')
    set(gca,'Xcolor','k','YColor','k','ZColor','k')
    view(-65,22);  axis on, hold on, pause(2),
    % boundIndex=boundIndex+1;  %Does not seem to be needed since FOR loop
    % automatically increases by 1 for every loop.
    clear EMG_SAMPLES EEG_SAMPLES EEG_TIMESTAMPS  %NOT IN ORGINAL
end

% Calculates the Mean and Std Deviation for sigma, theta and delta_lo values
%
[MeanValue,StdDevValue] = Calculate_MeanSTD_forPSDvalues(PSDvalues);
%
Sigma2sd_threshold = MeanValue.sigma + (2 * StdDevValue.sigma);
Sigma3sd_threshold = MeanValue.sigma + (3 * StdDevValue.sigma);
set(handles.sigma2SD_thresh,'String',num2str(Sigma2sd_threshold));
set(handles.sigma3SD_thresh,'String',num2str(Sigma3sd_threshold));

fid=fopen(newname,'a+');
fseek(fid,0,-1);  % Rewind the file completely and store these values
fprintf(fid,'The mean + 2*Std.dev value of Sigma power is:\t');
fprintf(fid,'%3.4f\t',Sigma2sd_threshold);
fprintf(fid,'The mean + 3*Std.dev value of Sigma power is:\t');
fprintf(fid,'%3.4f\n',Sigma3sd_threshold);
fclose(fid);

% Save the made figure to a fig file.
figure(1), figname1=strcat(pathname,'emg_st_dt_TR_',name1);
saveas(gcf,figname1,'fig');
figure(2),
saveas(gcf,figname,'fig');

figure(1),
X1=xlim; Y1=ylim; Z1=zlim; 
set(gca,'XlimMode','manual','YlimMode','manual','ZlimMode','manual')
figure(2),
X2=xlim; Y2=ylim; Z2=zlim;
set(gca,'XlimMode','manual','YlimMode','manual','ZlimMode','manual')

% Get the details in the Figure Manipulator GUI 
figureManipulations;
manipulator_handles=guihandles(figureManipulations);

set(manipulator_handles.x1mini,'String',X1(1)), set(manipulator_handles.x1maxi,'String',X1(2))
set(manipulator_handles.y1mini,'String',Y1(1)), set(manipulator_handles.y1maxi,'String',Y1(2))
set(manipulator_handles.z1mini,'String',Z1(1)), set(manipulator_handles.z1maxi,'String',Z1(2))
set(manipulator_handles.x2mini,'String',X2(1)), set(manipulator_handles.x2maxi,'String',X2(2))
set(manipulator_handles.y2mini,'String',Y2(1)), set(manipulator_handles.y2maxi,'String',Y2(2))
set(manipulator_handles.z2mini,'String',Z2(1)), set(manipulator_handles.z2maxi,'String',Z2(2))

fprintf(' Finished plotting the entire file... \n');

set(manipulator_handles.sig_del_theta_figure, 'String',strcat(figname,'.fig'));
set(manipulator_handles.emg_st_dt_figure, 'String',strcat(figname1,'.fig'));
guidata(manipulator_handles.FigureManipulations,manipulator_handles);

% --------------------------------------------------------------------
function autoscored_Callback(hObject, eventdata, handles)
global EMG_SAMPLES EEG_SAMPLES EEG_TIMESTAMPS EPOCHSIZE State INDEX EPOCHtime...
    Flag Statetime Statenum EPOCHnum Indices Cmap x1 y1 z1 

% This file uses PERSISTENT variables, and hence we clear them on every new use
clear('fft_psd_and_statescore_of_epoch')
%
filename1=get(handles.emgfile,'TooltipString');
filename2=get(handles.eegfile,'TooltipString');
filename3=get(handles.tstampsfile,'TooltipString');
filename4=get(handles.trainingfile,'TooltipString');
% Add to the search path of Matlab, the path where it can find the saved
% workspace variables. Remove the path at the end of the program.
if isempty(filename4) ==0 % If training filename is provided
    autoscoring_the_states_on_scatterplot(filename1,filename2,filename3,filename4);
else                    % carry on without the file
    autoscoring_the_states_on_scatterplot(filename1,filename2,filename3);
end
figure(1),xl=xlim;  yl=ylim ; zl = zlim;
set(handles.x1mini,'String',xl(1)), set(handles.x1maxi,'String',xl(2))
set(handles.y1mini,'String',yl(1)), set(handles.y1maxi,'String',yl(2))
set(handles.z1mini,'String',zl(1)), set(handles.z1maxi,'String',zl(2))
% --- Executes during object creation, after setting all properties.
function Sigma2SD_thresh_CreateFcn(hObject, eventdata, handles)
function Sigma2SD_thresh_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function sigma2SD_thresh_CreateFcn(hObject, eventdata, handles)
function sigma2SD_thresh_Callback(hObject, eventdata, handles)

%##########################################################################
%               CODE FOR SETTING DEFAULT VALUES OF FILTERS
%##########################################################################
% --- Executes on button press in Set_pushbutton.
function Set_pushbutton_Callback(hObject, eventdata, handles)
global D_lo D_hi T_lo T_hi S_lo S_hi B_lo B_hi
global EEG_Fc EMG_Fc EEG_highpass_enable EEG_Notch_enable EMG_Notch_enable EMG_lowpass_enable
global Input3_enable Input3_LP_enable Input3_Notch_enable Input3_HP_enable
global Input4_enable Input4_HP_enable Input4_Notch_enable Input4_LP_enable

D_lo = 0.4; set(handles.Delta_lo, 'String', D_lo);
D_hi = 4; set(handles.Delta_hi, 'String', D_hi);
T_lo = 5; set(handles.Theta_lo, 'String', T_lo);
T_hi = 9; set(handles.Theta_hi, 'String', T_hi);
S_lo = 10; set(handles.Sigma_lo, 'String', S_lo);
S_hi = 14; set(handles.Sigma_hi, 'String', S_hi);
B_lo = 15; set(handles.Beta_lo, 'String', B_lo);
B_hi = 20; set(handles.Beta_hi, 'String', B_hi);

EEG_Fc = 30; set(handles.EEG_cutoff, 'String', EEG_Fc);
EEG_highpass_enable=0; set(handles.EEG_Highpass_checkbox, 'Value', 0);
EEG_Notch_enable=0; set(handles.Notch60EEG, 'Value', 0);

EMG_Fc = 30; set(handles.EMG_cutoff, 'String', EMG_Fc);
EMG_Notch_enable=0; set(handles.Notch60EMG, 'Value', 0);
EMG_lowpass_enable=0; set(handles.EMG_LP_checkbox, 'Value', 0);

Input3_enable=0; set(handles.Input3_checkbox, 'Value', 0);
Input3_HP_enable=0; set(handles.Input3_HP_checkbox, 'Value', 0);
Input3_LP_enable=0; set(handles.Input3_LP_checkbox, 'Value', 0);
Input3_Notch_enable=0; set(handles.Notch60Input3, 'Value', 0);

Input4_enable=0; set(handles.Input4_checkbox, 'Value', 0);
Input4_HP_enable=0; set(handles.Input4_HP_checkbox, 'Value', 0);
Input4_LP_enable=0; set(handles.Input4_LP_checkbox, 'Value', 0);
Input4_Notch_enable=0; set(handles.Notch60Input4, 'Value', 0);

%##########################################################################
%           Code for User Input of EEG Wave Bandwidths
%##########################################################################
function Delta_lo_Callback(hObject, eventdata, handles)
global D_lo
Delta_lo = str2double(get(hObject,'String'));   % returns contents of Delta_lo as a double
if isnan(Delta_lo)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
D_lo = Delta_lo; 
% --- Executes during object creation, after setting all properties.
function Delta_lo_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Delta_hi_Callback(hObject, eventdata, handles)
global D_hi
Delta_hi = str2double(get(hObject,'String'));   % returns contents of EMG_Low as a double
if isnan(Delta_hi)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
D_hi = Delta_hi; 
% --- Executes during object creation, after setting all properties.
function Delta_hi_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Theta_lo_Callback(hObject, eventdata, handles)
global T_lo
Theta_lo = str2double(get(hObject,'String'));   % returns contents of Theta_lo as a double
if isnan(Theta_lo)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
T_lo = Theta_lo; 
% --- Executes during object creation, after setting all properties.
function Theta_lo_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Theta_hi_Callback(hObject, eventdata, handles)
global T_hi
Theta_hi = str2double(get(hObject,'String'));   % returns contents of Theta_hi as a double
if isnan(Theta_hi)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
T_hi = Theta_hi; 
% --- Executes during object creation, after setting all properties.
function Theta_hi_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Sigma_lo_Callback(hObject, eventdata, handles)
global S_lo
Sigma_lo = str2double(get(hObject,'String'));   % returns contents of Sigma_lo as a double
if isnan(Sigma_lo)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
S_lo = Sigma_lo;
% --- Executes during object creation, after setting all properties.
function Sigma_lo_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Sigma_hi_Callback(hObject, eventdata, handles)
global S_hi
Sigma_hi = str2double(get(hObject,'String'));   % returns contents of Sigma_hi as a double
if isnan(Sigma_hi)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
S_hi = Sigma_hi;
% --- Executes during object creation, after setting all properties.
function Sigma_hi_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Beta_lo_Callback(hObject, eventdata, handles)
global B_lo
Beta_lo = str2double(get(hObject,'String'));   % returns contents of Beta_lo as a double
if isnan(Beta_lo)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
B_lo = Beta_lo; 
% --- Executes during object creation, after setting all properties.
function Beta_lo_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Beta_hi_Callback(hObject, eventdata, handles)
global B_hi
Beta_hi = str2double(get(hObject,'String'));   % returns contents of Beta_hi as a double
if isnan(Beta_hi)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
B_hi = Beta_hi;
% --- Executes during object creation, after setting all properties.
function Beta_hi_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%##########################################################################
%                       User Input for EEG Filters
%##########################################################################
function EEG_cutoff_Callback(hObject, eventdata, handles)
global EEG_Fc
EEG_cutoff = str2double(get(hObject,'String'));   % returns contents of EEG_cutoff as a double
if isnan(EEG_cutoff)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if EEG_cutoff > 124
    %set(hObject, 'String', 0);
    errordlg('Cutoff frequency must be < 125 Hz','Error');
end
EEG_Fc = EEG_cutoff;
% 
function EEG_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in EEG_Highpass_checkbox_checkbox.
function EEG_Highpass_checkbox_Callback(hObject, eventdata, handles)
global EEG_highpass_enable
if (get(hObject,'Value') == get(hObject,'Max')) % then checkbox is checked-take approriate action
    EEG_highpass_enable=1;
else    % checkbox is not checked-take approriate action
    EEG_highpass_enable=0;
end

function EEG_Highpass_cutoff_Callback(hObject, eventdata, handles)
global EEG_HP_Fc
EEG_Highpass_cutoff = str2double(get(hObject,'String'));   % returns contents of EEG_cutoff as a double
if isnan(EEG_Highpass_cutoff)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if EEG_Highpass_cutoff < 0
    %set(hObject, 'String', 0);
    errordlg('Cutoff frequency must be a positive number','Error');
end
EEG_HP_Fc = EEG_Highpass_cutoff;

function EEG_Highpass_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in Notch60.
function Notch60EEG_Callback(hObject, eventdata, handles)
global EEG_Notch_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    EEG_Notch_enable=1;
else
    EEG_Notch_enable=0;
end

%##########################################################################
%                       User Input for EMG Filters
%##########################################################################
function EMG_cutoff_Callback(hObject, eventdata, handles)
global EMG_Fc
EMG_cutoff = str2double(get(hObject,'String'));   % returns contents of EMG_cutoff as a double
if isnan(EMG_cutoff)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
EMG_Fc = EMG_cutoff;
% --- Executes during object creation, after setting all properties.
function EMG_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in EMG_LP_checkbox_checkbox.
function EMG_LP_checkbox_Callback(hObject, eventdata, handles)
global EMG_lowpass_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    % then checkbox is checked-take approriate action
    EMG_lowpass_enable=1;
else
    % checkbox is not checked-take approriate action
    EMG_lowpass_enable=0;
end

function EMG_Lowpass_cutoff_Callback(hObject, eventdata, handles)
global EMG_LP_Fc
EMG_Lowpass_cutoff = str2double(get(hObject,'String'));   % returns contents of EMG_Lowpass_cutoff as a double
if isnan(EMG_Lowpass_cutoff)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if EMG_Lowpass_cutoff < 0
    %set(hObject, 'String', 0);
    errordlg('Cutoff frequency must be a positive number','Error');
end
EMG_LP_Fc = EMG_Lowpass_cutoff;

function EMG_Lowpass_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in Notch60EMG.
function Notch60EMG_Callback(hObject, eventdata, handles)
global EMG_Notch_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    % then checkbox is checked-take approriate action
    EMG_Notch_enable=1;
else
    % checkbox is not checked-take approriate action
    EMG_Notch_enable=0;
end

%##########################################################################
%                       OPTIONAL INPUT 3 Settings
%##########################################################################
function Input3file_Callback(hObject, eventdata, handles)
function Input3file_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in Input3_checkbox.
function Input3_checkbox_Callback(hObject, eventdata, handles)
global Input3_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    Input3_enable=1;
else
    Input3_enable=0;
end
% --- Executes on button press in Input3_HP_checkbox.
function Input3_HP_checkbox_Callback(hObject, eventdata, handles)
global Input3_HP_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    Input3_HP_enable=1;
else
    Input3_HP_enable=0;
end
%
function Input3_Highpass_cutoff_Callback(hObject, eventdata, handles)
global Input3_HP_Fc
Input3_Highpass_cutoff = str2double(get(hObject,'String'));   % returns contents of Input3_Highpass_cutoff as a double
if isnan(Input3_Highpass_cutoff)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if Input3_Highpass_cutoff < 0
    %set(hObject, 'String', 0);
    errordlg('Cutoff frequency must be a positive number','Error');
end
Input3_HP_Fc = Input3_Highpass_cutoff;
%
function Input3_Highpass_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in Input3_LP_checkbox.
function Input3_LP_checkbox_Callback(hObject, eventdata, handles)
global Input3_LP_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    Input3_LP_enable=1;
else
    Input3_LP_enable=0;
end
%
function Input3_Lowpass_cutoff_Callback(hObject, eventdata, handles)
global Input3_LP_Fc
Input3_Lowpass_cutoff = str2double(get(hObject,'String'));   % returns contents of Input3_Lowpass_cutoff as a double
if isnan(Input3_Lowpass_cutoff)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if Input3_Lowpass_cutoff < 0
    %set(hObject, 'String', 0);
    errordlg('Cutoff frequency must be a positive number','Error');
end
Input3_LP_Fc = Input3_Lowpass_cutoff;
% --- Executes during object creation, after setting all properties.
function Input3_Lowpass_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in Notch60Input3.
function Notch60Input3_Callback(hObject, eventdata, handles)
global Input3_Notch_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    Input3_Notch_enable=1;
else
    Input3_Notch_enable=0;
end
% --- Executes during object creation, after setting all properties.
function Notch60Input3_CreateFcn(hObject, eventdata, handles)

%##########################################################################
%                       OPTIONAL INPUT 4 Settings
%##########################################################################
function Input4file_Callback(hObject, eventdata, handles)
function Input4file_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in Input4_checkbox.
function Input4_checkbox_Callback(hObject, eventdata, handles)
global Input4_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    Input4_enable=1;
else
    Input4_enable=0;
end
% --- Executes on button press in Input4_HP_checkbox.
function Input4_HP_checkbox_Callback(hObject, eventdata, handles)
global Input4_HP_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    Input4_HP_enable=1;
else
    Input4_HP_enable=0;
end
%
function Input4_Highpass_cutoff_Callback(hObject, eventdata, handles)
global Input4_HP_Fc
Input4_Highpass_cutoff = str2double(get(hObject,'String'));   % returns contents of Input4_Highpass_cutoff as a double
if isnan(Input4_Highpass_cutoff)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if Input4_Highpass_cutoff < 0
    %set(hObject, 'String', 0);
    errordlg('Cutoff frequency must be a positive number','Error');
end
Input4_HP_Fc = Input4_Highpass_cutoff;
%
function Input4_Highpass_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in Input4_LP_checkbox.
function Input4_LP_checkbox_Callback(hObject, eventdata, handles)
global Input4_LP_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    % then checkbox is checked-take approriate action
    Input4_LP_enable=1;
else
    % checkbox is not checked-take approriate action
    Input4_LP_enable=0;
end
%
function Input4_Lowpass_cutoff_Callback(hObject, eventdata, handles)
global Input4_LP_Fc
Input4_Lowpass_cutoff = str2double(get(hObject,'String'));   % returns contents of Input4_Lowpass_cutoff as a double
if isnan(Input4_Lowpass_cutoff)
    %set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
if Input4_Lowpass_cutoff < 0
    %set(hObject, 'String', 0);
    errordlg('Cutoff frequency must be a positive number','Error');
end
Input4_LP_Fc = Input4_Lowpass_cutoff;
% --- Executes during object creation, after setting all properties.
function Input4_Lowpass_cutoff_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% --- Executes on button press in Notch60Input4.
function Notch60Input4_Callback(hObject, eventdata, handles)
global Input4_Notch_enable
if (get(hObject,'Value') == get(hObject,'Max'))
    Input4_Notch_enable=1;
else
    Input4_Notch_enable=0;
end
% --- Executes on selection change in SampfactorMenu.
function SampfactorMenu_Callback(hObject, eventdata, handles)
global sampfactor
val = get(hObject,'Value');
string_list = get(hObject,'String');
selected_string = string_list{val};
if strcmp(selected_string, 'Select')
    errordlg('Downsampling must be selected','Error');
else
    samp = str2double(selected_string);
end
sampfactor = samp;
% --- Executes during object creation, after setting all properties.
function SampfactorMenu_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tstampsfile_Callback(hObject, eventdata, handles)
% hObject    handle to tstampsfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tstampsfile as text
%        str2double(get(hObject,'String')) returns contents of tstampsfile as a double




% --------------------------------------------------------------------
function filemenu_Callback(hObject, eventdata, handles)
% hObject    handle to filemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
