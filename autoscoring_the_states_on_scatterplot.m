% This function decides the states based on the thresholds given in the main
% M-file of scorematic.m. It also generates the 3-D scatter plot of EMG,
% Sigma*Theta and Delta/Theta power levels per epoch

% This is linked with scorematic.m file
function[]=autoscoring_the_states_on_scatterplot(filename1,filename2,varargin)

global EMG_SAMPLES EEG_SAMPLES EEG_TIMESTAMPS EPOCHSIZE State INDEX Statetime...
    EPOCHtime Cmap Flag Statenum EPOCHnum X1 Y1 Z1 X2 Y2 Z2 stateTrack stateHistory
global exactLow exactHi boundIndex
handles=guihandles(scorematic);
% Take the variables out of the VARARGIN depending on its length.
if length(varargin) == 2
    filename3=char(varargin(1,1));
    filename4=char(varargin(1,2)); 
else
    filename3=char(varargin(1,1));
    filename4=[];
end

stateTrack = cell(1);
stateHistory = cell(1);
% Setup the colormap for the states.
Cmap(1,:)=[1 0.8 0];  % Yellow   => Active Waking
Cmap(2,:)=[0 0 1];    % Blue     => Quiet Sleep
Cmap(3,:)=[1 0 0];    % Red      => REM
Cmap(4,:)=[0 1 0.1];  % Green    => Quiet Waking
Cmap(5,:)=[0 0 0];    % Black    => Unhooked
Cmap(6,:)=[0 1 1];    % Cyan     => Transition to REM
Cmap(7,:)= [0.85 0.85 0.85];    % Clear    => Cleared State
Cmap(8,:)=[1 1 1];    % White    => Intermediate Waking

% Setup the colormap for the preview stage of a 3-D plot.
Pmap(1,:)=[1 0.5 0.5]; Pmap(2,:)=[0.6 0 0.8]; Pmap(3,:)=[1 0 1];
Pmap(4,:)=[.1 0.4 0.5];Pmap(5,:)=[0 0 0.3]; Pmap(6,:)=[0.2 0.5 1];
Pmap(7,:)=[0 1 0.5]; Pmap(8,:)=[0.8 0.2 0.1];
currentdate = date;
prompt={'Enter your name:','Enter the date you are scoring the file on:'};
def={'Brooks',currentdate};
dlgTitle='Input for file management';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
name=char(answer(1,:));
date1=char(answer(2,:));
% animal=char(answer(3,:));
% study=char(answer(4,:));

filename5=get(handles.autoscoredfile,'String');
SIGMAthresh=str2double(get(handles.sigma_thresh,'String'));
DELTAthresh=str2double(get(handles.delta_thresh,'String'));
THETAthresh=str2double(get(handles.theta_thresh,'String'));
STthresh=str2double(get(handles.st_thresh,'String'));
DTthresh=str2double(get(handles.dt_thresh,'String'));
EMGthresh=str2double(get(handles.emg_thresh,'String'));
Unhookedthresh=str2double(get(handles.unhooked_thresh,'String'));
Sigma3SDthresh=str2double(get(handles.sigma3SD_thresh,'String'));
Sigma2SDthresh=str2double(get(handles.sigma2SD_thresh,'String'));
working_dir=pwd;
current_dir='C:\SleepData\Results';
cd(current_dir);
[filename,pathname] = uiputfile(filename5,'Save Autoscored file as:');
filename=strcat(pathname,filename);
cd(working_dir);
% Make a strcuture of all thresholds
Thresholds = struct('sigma',SIGMAthresh,'delta',DELTAthresh,'theta',THETAthresh,...
    'sigmatheta',STthresh,'deltatheta',DTthresh,'emg',EMGthresh,...
    'unhooked',Unhookedthresh,'sigma3SD',Sigma3SDthresh,'sigma2SD',Sigma2SDthresh);

fid=fopen(filename,'w');
fprintf(fid,'Name:\t');              fprintf(fid,'%s\t',name);
fprintf(fid,'Date:\t');              fprintf(fid,'%s\t',date1);
fprintf(fid,'3*SD(sigmathresh):\t'); fprintf(fid,'%2.3d\t',Sigma3SDthresh);
fprintf(fid,'SigmaThresh:\t');       fprintf(fid,'%2.3d\n',SIGMAthresh);
% fprintf(fid,'Animal No:\t');         fprintf(fid,'%s\t',animal);
% fprintf(fid,'Study Info:\t');        fprintf(fid,'%s\t',study);
fprintf(fid,'DeltaThresh:\t');       fprintf(fid,'%2.3d\t',DELTAthresh);
fprintf(fid,'ThetaThresh:\t');       fprintf(fid,'%2.3d\n',THETAthresh);
fprintf(fid,'ST Thresh:\t');         fprintf(fid,'%2.3d\t',STthresh);
fprintf(fid,'DT Thresh:\t');         fprintf(fid,'%2.3d\t',DTthresh);
fprintf(fid,'EMGthresh:\t');         fprintf(fid,'%2.3d\t',EMGthresh);
fprintf(fid,'Unhookedthresh:\t');    fprintf(fid,'%2.3d\n',Unhookedthresh);
fprintf(fid,'Index');                fprintf(fid,'\t');         
fprintf(fid,'Timestamp ');           fprintf(fid,'\t');
fprintf(fid,'State\n');
fclose(fid);

[tbounds]=xlsread(filename3);
lbound=tbounds(1:end,1);
ubound=tbounds(1:end,2);
exactLowIndx = tbounds(1:end,3); 
exactHiIndx = tbounds(1:end,4);
Flag=1;  % Indicates that fft_calc_and_statescore_of_epoch is used for auto scoring purposes.
INDEX=1;  Index_start=1;
figure(1);     whitebg(figure(1),[.165,.465,.665]), set(gcf,'position',[672 533 560 420])
figure(2);     whitebg(figure(2),[.165,.465,.665]), set(gcf,'position',[672 33 560 420])

for boundIndex=1:length(ubound)
    lowerbound=lbound(boundIndex);
    upperbound=ubound(boundIndex);
    exactLow = exactLowIndx(boundIndex);
    exactHi = exactHiIndx(boundIndex);
    %fprintf('          Starting to process the iteration #%d ......\n',nBoundIndices);
    
    if isempty(filename4) ==0 % If training filename is provided
        [PSDvalues] = read_n_extract_datafiles(Thresholds,filename1,filename2,filename3,filename4,num2str(lowerbound),num2str(upperbound));
    else                    % carry on without the file
        [PSDvalues] = read_n_extract_datafiles(Thresholds,filename1,filename2,filename3,num2str(lowerbound),num2str(upperbound));
    end
    fprintf(' Finished fileread .... \n');
    S=[];C=[]; awind=[];qsind=[];reind=[];qwind=[];unind=[];iwind=[];
    S=ones(length(PSDvalues.deltatheta),1);  S(:)=12;
    C=zeros(length(PSDvalues.deltatheta),3);
    
    awind=find(Statenum == 1);
    C(awind,1)=Cmap(1,1);C(awind,2)=Cmap(1,2);C(awind,3)=Cmap(1,3);
    qsind=find(Statenum == 2);
    C(qsind,1)=Cmap(2,1);C(qsind,2)=Cmap(2,2);C(qsind,3)=Cmap(2,3);
    reind=find(Statenum == 3);
    C(reind,1)=Cmap(3,1);C(reind,2)=Cmap(3,2);C(reind,3)=Cmap(3,3);
    qwind=find(Statenum == 4);
    C(qwind,1)=Cmap(4,1);C(qwind,2)=Cmap(4,2);C(qwind,3)=Cmap(4,3);
    unind=find(Statenum == 5);
    C(unind,1)=Cmap(5,1);C(unind,2)=Cmap(5,2);C(unind,3)=Cmap(5,3);
    trind=find(Statenum == 6);
    C(trind,1)=Cmap(6,1);C(trind,2)=Cmap(6,2);C(trind,3)=Cmap(6,3);
    iwind=find(Statenum == 8);
    C(iwind,1)=Cmap(8,1);C(iwind,2)=Cmap(8,2);C(iwind,3)=Cmap(8,3);
    %fprintf(' Epoch power related with thresholds....\n');  
    figure(1),scatter3(PSDvalues.deltatheta,PSDvalues.sigmatheta,...
        PSDvalues.emg,S,C,'filled');
    xlabel('Delta/Theta Power');ylabel('Sigma*Theta Power');zlabel('EMG Power');
    view(-65,22); axis on, hold on, pause(2),
    
    figure(2),scatter3(PSDvalues.sigma,PSDvalues.theta,PSDvalues.delta,S,C,'filled');
    xlabel('Sigma Power');ylabel('Theta Power');zlabel('Delta Power');
    axis on, hold on, pause(2),view(-65,22);
        
    Index_track=Index_start:Index_start+(INDEX-2);
    i=1;
    fid=fopen(filename,'a+');
    while i <= length(State)
        fprintf(fid,'%d',Index_track(i));   % The corresponding INDEX
        fprintf(fid,'\t');
        fprintf(fid,'%3.4f',Statetime(i));  % Time Stamp values
        fprintf(fid,'\t');
        fprintf(fid,'%6s', State(i,:));     % State scored automatically
        fprintf(fid,'\t');
        fprintf(fid,'%s', stateTrack{i});   % State tracker final state
        fprintf(fid,'\t');
        fprintf(fid,'%s', stateHistory{i});   % State tracker history
        fprintf(fid,'\n');
        i=i+1;
    end
    fclose(fid);
    Index_start=Index_start+INDEX-1;    
end    

figure(1),X1=xlim; Y1=ylim; Z1=zlim; xlim('manual'),ylim('manual'),zlim('manual');
figure(2),X2=xlim; Y2=ylim; Z2=zlim; xlim('manual'),ylim('manual'),zlim('manual');

% Save the made figure to a fig file.

[path,name,ext]=fileparts(filename5);
figure(1),  figname1=strcat(pathname,'emg_dt_st_',name);
saveas(gcf,figname1,'fig');
figure(2),  figname2=strcat(pathname,'sdt_',name);
saveas(gcf,figname2,'fig');

% Get the details in the Figure Manipulator GUI 
figureManipulations
manipulator_handles=guihandles(figureManipulations);
set(manipulator_handles.x1mini,'String',X1(1)), set(manipulator_handles.x1maxi,'String',X1(2))
set(manipulator_handles.y1mini,'String',Y1(1)), set(manipulator_handles.y1maxi,'String',Y1(2))
set(manipulator_handles.z1mini,'String',Z1(1)), set(manipulator_handles.z1maxi,'String',Z1(2))
set(manipulator_handles.x2mini,'String',X2(1)), set(manipulator_handles.x2maxi,'String',X2(2))
set(manipulator_handles.y2mini,'String',Y2(1)), set(manipulator_handles.y2maxi,'String',Y2(2))
set(manipulator_handles.z2mini,'String',Z2(1)), set(manipulator_handles.z2maxi,'String',Z2(2))

fprintf(' Finished scoring the entire file... \n');
set(manipulator_handles.emg_st_dt_figure,...
    'String',strcat(figname1,'.fig'))
set(manipulator_handles.sig_del_theta_figure,...
    'String',strcat(figname2,'.fig'))
guidata(manipulator_handles.FigureManipulations,manipulator_handles);
