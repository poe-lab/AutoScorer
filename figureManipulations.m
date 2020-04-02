function varargout = figureManipulations(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @figureManipulations_OpeningFcn, ...
                   'gui_OutputFcn',  @figureManipulations_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before figureManipulations is made visible.
function figureManipulations_OpeningFcn(hObject, eventdata, figurehandles, varargin)
% Choose default command line output for figureManipulations
figurehandles.output = hObject;
% Update figurehandles structure
guidata(hObject, figurehandles);

% --- Outputs from this function are returned to the command line.
function varargout = figureManipulations_OutputFcn(hObject, eventdata, figurehandles)
% Get default command line output from figurehandles structure
varargout{1} = figurehandles.output;

% --- Executes on button press in emgDeltathetaView.
function emgDeltathetaView_Callback(hObject, eventdata, figurehandles)
global X1 Y1 Z1 X2 Y2 Z2

figure(2), xlim([X2(1) X2(2)]), 
ylim([Y2(1) Y2(2)]), zlim([Z2(1) Z2(2)]), view(0,0)
figure(1), xlim([X1(1) X1(2)]), 
ylim([Y1(1) Y1(2)]), zlim([Z1(1) Z1(2)]), view( 0,0)

% --- Executes on button press in SigthetaDelthetaView.
function SigthetaDelthetaView_Callback(hObject, eventdata, figurehandles)
global X1 Y1 Z1 X2 Y2 Z2

figure(1), xlim([X1(1) X1(2)]), 
ylim([Y1(1) Y1(2)]), zlim([Z1(1) Z1(2)]), view( 0,90)
figure(2), xlim([X2(1) X2(2)]), 
ylim([Y2(1) Y2(2)]), zlim([Z2(1) Z2(2)]), view( 0,90)

% --- Executes on button press in emg_sigmathetaview.
function emg_sigmathetaview_Callback(hObject, eventdata, handles)
global X1 Y1 Z1 X2 Y2 Z2

figure(1), xlim([X1(1) X1(2)]), 
ylim([Y1(1) Y1(2)]), zlim([Z1(1) Z1(2)]), view( -90,0)
figure(2), xlim([X2(1) X2(2)]), 
ylim([Y2(1) Y2(2)]), zlim([Z2(1) Z2(2)]), view( -90,0)

% --- Executes on button press in restoreOriginalView.
function restoreOriginalView_Callback(hObject, eventdata, figurehandles)
global X1 Y1 Z1 X2 Y2 Z2 

figure(1), xlim([X1(1) X1(2)]), 
ylim([Y1(1) Y1(2)]), zlim([Z1(1) Z1(2)]), view( -65,22)
set(figurehandles.x1mini,'String',X1(1)), set(figurehandles.x1maxi,'String',X1(2))
set(figurehandles.y1mini,'String',Y1(1)), set(figurehandles.y1maxi,'String',Y1(2))
set(figurehandles.z1mini,'String',Z1(1)), set(figurehandles.z1maxi,'String',Z1(2))

figure(2), xlim([X2(1) X2(2)]), 
ylim([Y2(1) Y2(2)]), zlim([Z2(1) Z2(2)]), view( -65,22)
set(figurehandles.x2mini,'String',X2(1)), set(figurehandles.x2maxi,'String',X2(2))
set(figurehandles.y2mini,'String',Y2(1)), set(figurehandles.y2maxi,'String',Y2(2))
set(figurehandles.z2mini,'String',Z2(1)), set(figurehandles.z2maxi,'String',Z2(2))

% --- Executes on button press in refreshLimitsView.
function refreshLimitsView_Callback(hObject, eventdata, figurehandles)

x1(1)=str2double(get(figurehandles.x1mini,'string')); x1(2)=str2double(get(figurehandles.x1maxi,'string'));
y1(1)=str2double(get(figurehandles.y1mini,'string')); y1(2)=str2double(get(figurehandles.y1maxi,'string'));
z1(1)=str2double(get(figurehandles.z1mini,'string')); z1(2)=str2double(get(figurehandles.z1maxi,'string'));
figure(1), xlim([x1(1) x1(2)]), ylim([y1(1) y1(2)]), zlim([z1(1) z1(2)])

x2(1)=str2double(get(figurehandles.x2mini,'string')); x2(2)=str2double(get(figurehandles.x2maxi,'string'));
y2(1)=str2double(get(figurehandles.y2mini,'string')); y2(2)=str2double(get(figurehandles.y2maxi,'string'));
z2(1)=str2double(get(figurehandles.z2mini,'string')); z2(2)=str2double(get(figurehandles.z2maxi,'string'));
figure(2), xlim([x2(1) x2(2)]), ylim([y2(1) y2(2)]), zlim([z2(1) z2(2)])

% --------------------------------------------------------------------
function emg_st_dt_open_Callback(hObject, eventdata, figurehandles)
working_dir=pwd;
current_dir='C:\1-SAT2\Results\';
cd(current_dir);
[filename, pathname] = uigetfile('*.fig', 'Pick the EMG-ST-DT figure ');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    emg_st_dt_file= fullfile(pathname, filename);
    set(figurehandles.emg_st_dt_figure,'string',filename);
    set(figurehandles.emg_st_dt_figure,'Tooltipstring',emg_st_dt_file);
end

% --------------------------------------------------------------------
function sig_del_theta_open_Callback(hObject, eventdata, figurehandles)
working_dir=pwd;
current_dir='C:\1-SAT2\Results\';
cd(current_dir);
[filename, pathname] = uigetfile('*.fig', 'Pick the Sig-Del_Theta figure ');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    sig_del_theta_file= fullfile(pathname, filename);
    set(figurehandles.sig_del_theta_figure,'string',filename);
    set(figurehandles.sig_del_theta_figure,'Tooltipstring',sig_del_theta_file);
end

% --------------------------------------------------------------------
function figures_load_Callback(hObject, eventdata, figurehandles)
global X1 Y1 Z1 X2 Y2 Z2 
% Code for loading and opening the figures
fig1=get(figurehandles.emg_st_dt_figure,'Tooltipstring');
fig2=get(figurehandles.sig_del_theta_figure,'Tooltipstring');
if isempty(fig1) == 0
    open(fig1);
    figure(1),
    X1=xlim; Y1=ylim ; Z1=zlim; xlim('manual'),ylim('manual'),zlim('manual')
    set(figurehandles.x1mini,'String',X1(1)), set(figurehandles.x1maxi,'String',X1(2))
    set(figurehandles.y1mini,'String',Y1(1)), set(figurehandles.y1maxi,'String',Y1(2))
    set(figurehandles.z1mini,'String',Z1(1)), set(figurehandles.z1maxi,'String',Z1(2))

end
if isempty(fig2) == 0
    open(fig2);
    figure(2),
    X2=xlim;Y2=ylim;Z2=zlim;  xlim('manual'),ylim('manual'),zlim('manual');
    set(figurehandles.x2mini,'String',X2(1)), set(figurehandles.x2maxi,'String',X2(2))
    set(figurehandles.y2mini,'String',Y2(1)), set(figurehandles.y2maxi,'String',Y2(2))
    set(figurehandles.z2mini,'String',Z2(1)), set(figurehandles.z2maxi,'String',Z2(2))
end

% --------------------------------------------------------------------
function psdvaluesfile_open_Callback(hObject, eventdata, figurehandles)
global MeanValue StdDevValue

sigma=[];theta=[];delta=[];%npoints=[];
working_dir=pwd;
current_dir='C:\1-SAT2\Results\';
cd(current_dir);
[filename, pathname] = uigetfile('*.xls', 'Pick the PSDvalues(FFT) file');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please press the button again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    psdfile= fullfile(pathname, filename);
    copy_expression=['copy ' psdfile ' ' working_dir '\' filename ]; 
    dos(copy_expression);
    % Open and read the PSD values from the FFTfile just selected 
    [psdvalues ,columnname] = xlsread_modified(filename);
    del_expression=['del ' working_dir '\' filename ];
    dos(del_expression);
    
    sigma = psdvalues(1:end-2,2);
    theta = psdvalues(1:end-2,3);
    delta = psdvalues(1:end-2,4);
    
    sumsigma=sum(sigma);
    sumtheta=sum(theta);
    sumdelta=sum(delta);
    squaresumsigma=sum(sigma.^2);
    squaresumtheta=sum(theta.^2);
    squaresumdelta=sum(delta.^2);
    nPointssigma = length(sigma);
    nPointstheta = length(theta);
    nPointsdelta = length(delta);
    
    squaresum=struct('delta',squaresumdelta, 'theta',squaresumtheta,...
        'sigma',squaresumsigma);
    SuM=struct('delta',sumdelta, 'theta',sumtheta, 'sigma',sumsigma);
    nPoints=struct('delta',nPointsdelta, 'theta',nPointstheta,...
        'sigma',nPointssigma);
    PSDvalues = struct( 'delta',delta, 'theta',theta,'sigma',sigma, ...
        'squaresum',squaresum,'sum',SuM, 'nPoints',nPoints);
    
    % Calculates the Mean and Std Deviation for sigma, theta and delta values
    %
    [MeanValue,StdDevValue] = Calculate_MeanSTD_forPSDvalues(PSDvalues);
    %
end

% --------------------------------------------------------------------- 
function SD1view_Callback(hObject, eventdata, figurehandles)
global MeanValue StdDevValue X1 Y1 Z1 X2 Y2 Z2
persistent I1 I2 I3 J1 J2 J3 K1 K2 K3 L1 L2 L3
if get(figurehandles.SD1view,'value') == 1
    % Plot the Std Dev. values on figure(2)
    
    xref = 0:X2(2); yref = 0:Y2(2); zref = 0:Z2(2);
    array=[X2(2) Y2(2) Z2(2)]; [maxarray,I]=max(array);
    zeroref=zeros(1,array(I));
    figure(2),
    I1=plot3(repmat(StdDevValue.sigma,1,Y2(2)),yref(2:end),zeroref(1:Y2(2)),'y');
    J1=plot3(repmat(StdDevValue.sigma,1,Y2(2)),yref(2:end),repmat(Z2(2),1,Y2(2)),'y');
    K1=plot3(repmat(StdDevValue.sigma,1,Z2(2)),zeroref(1:Z2(2)),zref(2:end),'y');
    L1=plot3(repmat(StdDevValue.sigma,1,Z2(2)),repmat(Y2(2),1,Z2(2)),zref(2:end),'y');
    
    I2=plot3(xref(2:end),repmat(StdDevValue.theta,1,X2(2)),zeroref(1:X2(2)),'k');
    J2=plot3(xref(2:end),repmat(StdDevValue.theta,1,X2(2)),repmat(Z2(2),1,X2(2)),'k');
    K2=plot3(zeroref(1:Z2(2)),repmat(StdDevValue.theta,1,Z2(2)),zref(2:end),'k');
    L2=plot3(repmat(X2(2),1,Z2(2)),repmat(StdDevValue.theta,1,Z2(2)),zref(2:end),'k');
    
    I3=plot3(zeroref(1:Y2(2)),yref(2:end),repmat(StdDevValue.delta,1,Y2(2)),'c');
    J3=plot3(repmat(X2(2),1,Y2(2)),yref(2:end),repmat(StdDevValue.delta,1,Y2(2)),'c');
    K3=plot3(xref(2:end),zeroref(1:X2(2)),repmat(StdDevValue.delta,1,X2(2)),'c');
    L3=plot3(xref(2:end),repmat(Y2(2),1,X2(2)),repmat(StdDevValue.delta,1,X2(2)),'c');
    
    legend([L1,L2,L3],'Sigma','Theta','Delta',1);
else
    set([I1,J1,K1,L1],'visible','off');
    set([I2,J2,K2,L2],'visible','off');
    set([I3,J3,K3,L3],'visible','off');
end

% --- Executes on button press in SD2view.
function SD2view_Callback(hObject, eventdata, figurehandles)
% --- Executes on button press in SD3view.
function SD3view_Callback(hObject, eventdata, figurehandles)
% --------------------------------------------------------------------
function filemenu_Callback(hObject, eventdata, figurehandles)
% --------------------------------------------------------------------
function load_menu_Callback(hObject, eventdata, figurehandles)

