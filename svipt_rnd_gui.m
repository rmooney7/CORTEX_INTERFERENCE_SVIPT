function varargout = salas_gui91(varargin)
% AXESMENUTOOLBAR Example for creating GUIs with menus, toolbar, and plots
%       AXESMENUTOOLBAR is an example GUI for demonstrating how to creating
%       GUIs using nested functions. It shows how to generate different
%       plots and how to add menus and toolbar to the GUIs.

%   Copyright 1984-2006 The MathWorks, Inc.

% Declare non-UI data so that they can be used in any functions in this GUI
% file, including functions triggered by creating the GUI layout below

yScale=30/32768; %for use with the CED, conversion to voltage

mOutputArgs     = {};       % Variable for storing output when GUI returns

%define the types of functions that can appear in the function menu
mFunctionTypes      = {'LI','E','LO','R','CI'};

% Declare and create all the UI objects in this GUI here so that they can
% be used in any functions
hMainFigure    =   figure(...       % the main GUI figure
    'MenuBar','none', ...
    'Toolbar','none', ...
    'HandleVisibility','callback', ...
    'Name', mfilename, ...
    'NumberTitle','off', ...
    'Color', get(0, 'defaultuicontrolbackgroundcolor'));
hPlotAxes      =   axes(...         % the axes for plotting selected plot
    'Parent', hMainFigure, ...
    'Units', 'normalized', ...
    'HandleVisibility','callback', ...
    'Position',[0.10 0.10 0.8 0.4]);
hTrialDoneAxes      =   axes(...         % the axes for plotting selected plot
    'Parent', hMainFigure, ...
    'Units', 'normalized', ...
    'HandleVisibility','callback', ...
    'Position',[0.45 0.15 0.1 0.1]);
hTextAxes      =   axes(...         % the axes for plotting selected plot
    'Parent', hMainFigure, ...
    'Units', 'normalized', ...
    'HandleVisibility','callback', ...
    'Position',[0.10 0.4 0.8 0.1]);
hReadyAxes      =   axes(...         % the axes for plotting selected plot
    'Parent', hMainFigure, ...
    'Units', 'normalized', ...
    'HandleVisibility','callback', ...
    'Position',[0.75 0.60 0.15 0.2]);
hFunctionPopupmenu =   uicontrol(...    % list of available types of plot
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'Position',[0.3 0.7 0.1 0.05],...
    'HandleVisibility','callback', ...
    'String',mFunctionTypes,...
    'Style','popupmenu');
hFunctionPopupmenuText  =   uicontrol(...    % box that holds where to save files
    'Style','text',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.3 0.75 0.5 0.05],...
    'String','Function Type',...
    'HorizontalAlignment','left');
hStartButton  =   uicontrol(...    % button for updating selected plot
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.7 0.85 0.1 0.05],...
    'String','START',...
    'Callback', @hStartButtonCallback);
hTrialStartButton  =   uicontrol(...    % button for updating selected plot
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.8 0.85 0.1 0.05],...
    'String','Pause',...
    'Visible','on',...
    'Callback', @hTrialStartButtonCallback);
hFilebox  =   uicontrol(...    % box that holds where to save files
    'Style','edit',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.10 0.90 0.5 0.05],...
    'String','Enter Filename',...
    'HorizontalAlignment','left');
hFileboxText  =   uicontrol(...    % box that holds where to save files
    'Style','text',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.10 0.95 0.5 0.05],...
    'String','Enter Output Filename',...
    'HorizontalAlignment','left');
hGain1box  =   uicontrol(...    % box that holds where to save files
    'Style','edit',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.45 0.7 0.05 0.05],...
    'String','1',...
    'HorizontalAlignment','left');
hGain1boxText  =   uicontrol(...    % box that holds where to save files
    'Style','text',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.45 0.75 0.05 0.05],...
    'String','a',...
    'HorizontalAlignment','left');
hGain2box  =   uicontrol(...    % box that holds where to save files
    'Style','edit',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.55 0.7 0.05 0.05],...
    'String','1',...
    'HorizontalAlignment','left');
hGain2boxText  =   uicontrol(...    % box that holds where to save files
    'Style','text',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.55 0.75 0.05 0.05],...
    'String','b',...
    'HorizontalAlignment','left');
hGain3box  =   uicontrol(...    % box that holds where to save files
    'Style','edit',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.65 0.7 0.05 0.05],...
    'String','0',...
    'HorizontalAlignment','left');
hGain3boxText  =   uicontrol(...    % box that holds where to save files
    'Style','text',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.65 0.75 0.05 0.05],...
    'String','c',...
    'HorizontalAlignment','left');
hTimerbox  =   uicontrol(...    % box that holds where to save files
    'Style','edit',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.10 0.8 0.15 0.05],...
    'String','0',...
    'HorizontalAlignment','left');
hTimerText  =   uicontrol(...    % box that holds where to save files
    'Style','text',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.10 0.85 0.15 0.05],...
    'String','Timer',...
    'HorizontalAlignment','left');
hTrialbox  =   uicontrol(...    % box that holds where to save files
    'Style','edit',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.3 0.8 0.05 0.05],...
    'String','0',...
    'HorizontalAlignment','left');
hTrialText  =   uicontrol(...    % box that holds where to save files
    'Style','text',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.30 0.85 0.15 0.05],...
    'String','# Trials',...
    'HorizontalAlignment','left');
hLagbox  =   uicontrol(...    % box that holds where to save files
    'Style','edit',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.5 0.8 0.15 0.05],...
    'String','0',...
    'HorizontalAlignment','left');
hLagText  =   uicontrol(...    % box that holds where to save files
    'Style','text',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.50 0.85 0.15 0.05],...
    'String','% Lag',...
    'HorizontalAlignment','left');
hMRIcheck  =   uicontrol(...    % box that holds where to save files
    'Style','checkbox',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.10 0.75 0.2 0.05],...
    'String','Trigger MRI',...
    'Callback', @hMRIcheckCallback);
hTMScheck  =   uicontrol(...    % box that holds where to save files
    'Style','checkbox',...
    'Parent', hMainFigure, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.10 0.7 0.2 0.05],...
    'String','Trigger TMS',...
    'Callback', @hTMScheckCallback);
hFileMenu      =   uimenu(...       % File menu
    'Parent',hMainFigure,...
    'HandleVisibility','callback', ...
    'Label','File');
hOpenMenuitem  =   uimenu(...       % Open menu item
    'Parent',hFileMenu,...
    'Label','Open',...
    'HandleVisibility','callback', ...
    'Callback', @hOpenMenuitemCallback);
hPrintMenuitem  =  uimenu(...       % Print menu item
    'Parent',hFileMenu,...
    'Label','Print',...
    'HandleVisibility','callback', ...
    'Callback', @hPrintMenuitemCallback);
hCloseMenuitem  =  uimenu(...       % Close menu item
    'Parent',hFileMenu,...
    'Label','Close',...
    'Separator','on',...
    'HandleVisibility','callback', ...
    'Callback', @hCloseMenuitemCallback);
hToolbar       =   uitoolbar(...    % Toolbar for Open and Print buttons
    'Parent',hMainFigure, ...
    'HandleVisibility','callback');
hOpenPushtool  =   uipushtool(...   % Open toolbar button
    'Parent',hToolbar,...
    'TooltipString','Open File',...
    'CData',iconRead(fullfile(matlabroot, '/toolbox/matlab/icons/opendoc.mat')),...
    'HandleVisibility','callback', ...
    'ClickedCallback', @hOpenMenuitemCallback);
hPrintPushtool =   uipushtool(...    % Print toolbar button
    'Parent',hToolbar,...
    'TooltipString','Print Figure',...
    'CData',iconRead(fullfile(matlabroot, '/toolbox/matlab/icons/printdoc.mat')),...
    'HandleVisibility','callback', ...
    'ClickedCallback', @hPrintMenuitemCallback);


%get all initialization paremeters from the .ini file
params.ini_Filename = strcat(mfilename,'.ini');
fid = fopen(params.ini_Filename);
while(true)
    nameline = fgetl(fid);
    if(nameline == -1)
        break;
    elseif(strcmp(nameline, 'default_Filename'))
        valueline = fgetl(fid);
        params.(nameline) = valueline;
    else
        valueline = fgetl(fid);
        params.(nameline) = str2num(valueline);
    end
end
params;
fclose(fid);

%set all the required properties
set(hFilebox,'String',params.default_Filename);
set(hTrialbox,'String',params.default_NTrials);


%set all the current statuses
global StartButton_status;
StartButton_status = 0;
global TMStrigger_status;
TMStrigger_status = 0;
global MRItrigger_status;
MRItrigger_status = 0;
global TrialStartButton_status;
TrialStartButton_status = 1;

global ntrials;
global selected_function;
global a;
global b;
global c;
global tms_trigger;
global mri_trigger;

mouseover = 0;

tp=0;

x=1:100;
beep1=sin(1000*x);
beep2=sin(7500*x);

%this sets up the initial ready buddton
hbutt = rectangle('Parent',hReadyAxes,'Position',[-0.5,-0.5,1,1],'FaceColor','r','Curvature',1);
axis(hReadyAxes,[-1,1,-1,1],'off');

axis(hTextAxes,[0,1,0,1],'off');
for i = 1:size(params.text_locs,2)
    text(params.text_locs(i),0.5,int2str(params.text_data(i)),'Parent',hTextAxes);
end

%this sets up the initial ready buddton
hTrialDone = rectangle('Parent',hTrialDoneAxes,'Position',[-0.25,-0.25,0.5,0.5],'FaceColor','r','Curvature',1,'LineStyle','none');
axis(hTrialDoneAxes,[-1,1,-1,1],'off');

%%%%%

% This sets up the initial plot with the bars
for i = 1:size(params.checkpoint_locs,2)
    hCheckpoint(i) = rectangle('Parent',hPlotAxes,'Position',[params.checkpoint_locs(i)-params.checkpoint_thickness(i)/2,params.checkpoint_center_loc(i)-params.checkpoint_height(i)/2,params.checkpoint_thickness(i),params.checkpoint_height(i)],'FaceColor',[params.checkpoint_color_red(i),params.checkpoint_color_green(i),params.checkpoint_color_blue(i)],'LineStyle','none');
end

% % This sets up the initial plot with the bars AS BOXES (facecolor is set to 0; linecolor is set to original facecolor)
% for i = 1:size(params.checkpoint_locs,2)
%     hCheckpoint(i) = rectangle('Parent',hPlotAxes,'Position',[params.checkpoint_locs(i)-params.checkpoint_thickness(i)/2,params.checkpoint_center_loc(i)-params.checkpoint_height(i)/2,params.checkpoint_thickness(i),params.checkpoint_height(i)],'EdgeColor',[params.checkpoint_color_red(i),params.checkpoint_color_green(i),params.checkpoint_color_blue(i)],'LineStyle','-','LineWidth',2);
% end

%%%%

%setup the square dot
hRect = rectangle('Parent',hPlotAxes,'Position',[0-params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)],'FaceColor',params.pinchcursor_color);
axis(hPlotAxes,[0,1,0,1]);


%display axes - used in calibrating
if (params.display_axes==0)
    axis(hPlotAxes,'off');
else
    axis(hPlotAxes,'on');
end

%set the default function type
set(hGain1box,'String',num2str(params.a_init));
set(hGain2box,'String',num2str(params.b_init));
set(hGain3box,'String',num2str(params.c_init));
set(hFunctionPopupmenu,'Value',params.function_type);





%open the DAC Box
% yScale=25/32768; % results in Volts
res=matced64c('cedOpenX',0);  % open a connection
if (res < 0)
    disp(['1401 not opened, error number ' int2str(res)]);
else %opened OK
    res=matced64c('cedResetX'); %reset it
end
matced64c('cedLdX','C:\1401\','ADCMEM','ADCBST');

res=matced64c('cedWorkingSet',400,4000);
if (res > 0)
    disp('error with call to cedWorkingSet - try commenting it out');
    return
end


% Define default output and return it if it is requested by users
mOutputArgs{1} = hMainFigure;
if nargout>0
    [varargout{1:nargout}] = mOutputArgs{:};
end

%----------------------------------------------------------------------
    function hStartButtonCallback(hObject, eventdata)
        % Callback function run when the update button is pressed


        if(StartButton_status == 0)

            %change signal
            set(hbutt,'FaceColor','g');

            %open file
            filename = get(hFilebox,'String');
            fid = fopen(filename,'wt+');
            
            %create file that contains raw data
            filename_raw = strcat(filename,'raw');
            fid_raw = fopen(filename_raw,'wt+');

            %start experiment
            %get the function properties
            ntrials = str2num(get(hTrialbox,'String'));
            selected_function = get(hFunctionPopupmenu,'value');
            a = str2num(get(hGain1box,'String'));
            b = str2num(get(hGain2box,'String'));
            c = str2num(get(hGain3box,'String'));
            tms_trigger = get(hTMScheck,'value');
            mri_trigger = get(hMRIcheck,'value');

            %make the boxes invisible
            set(hFunctionPopupmenu,'Visible','off');
            set(hGain1box,'Visible','off');
            set(hGain2box,'Visible','off');
            set(hGain3box,'Visible','off');
            set(hLagbox,'Visible','off');
            set(hTMScheck,'Visible','off');
            set(hMRIcheck,'Visible','off');
            set(hFunctionPopupmenuText,'Visible','off');
            set(hGain1boxText,'Visible','off');
            set(hGain2boxText,'Visible','off');
            set(hGain3boxText,'Visible','off');
            set(hLagText,'Visible','off');

            %change the start button status
            StartButton_status;
            set(hStartButton,'String','STOP');
            StartButton_status = 1;

            %take down initial time
            tinit = clock;

            random_function         = load('random.txt');
            %interleaved_function    = load('interleaved.txt');

            %begin trials
            tp=0;
            tic
            for trialcount = 1:ntrials

                %set status of beginning of trial
                TrialStartButton_status = 1;



                %set the current status of the current location of the
                %cursor
                cursor_loc = 1;

                %set the flag that tells whether we are within the
                %checkpoints
                stays_in = 0;

                %display in trial box
                set(hTrialbox,'String',int2str(trialcount));
                
                %turn lag box off
                set(hLagbox,'Visible','off');
                set(hLagText,'Visible','off');
                
                %reset trial done circle
                set(hTrialDone,'FaceColor','g');

                %initialize count vector
                pinch_data_count = 0;
                %measure of lag
                lag_count = 0;

               
                tstart = (toc-tp);
                
                %main data acquisition and recording loop
                while StartButton_status == 1




                   

                    % % %
                    %better dampling of data
                    %matced64c('cedLdX','C:\1401\','ADCMEM','ADCBST');
                    res = -1;
                    matced64c('cedSendString','ADCBST,I,2,0,20,0,1,H,32,25;');
                    while res ~=0
                        matced64c('cedSendString','ADCBST,?;');
                        res=eval(matced64c('cedGetString'));
                    end
                    
                    pinch_data_temp1=mean(matced64c('cedToHost',10,0))*yScale;
                    
                    %Mutliply voltage by 2 when using new IPM650 to account
                    %for 0-5 V output instead of 0-10 V
%                     pinch_data_temp1=mean(matced64c('cedToHost',10,0))*yScale*2;

%Filter out valutes under .05 (eliminates visible baseline noise)
                    if pinch_data_temp1<.05
                        pinch_data_temp1=0;
                    end
                    pinch_data_temp=pinch_data_temp1;



                    %generate fake data
%                     pinch_data_temp=sin(pinch_data_count/50);
%                     timeper = 10;
%                     if ceil(pinch_Fdata_count/(100*timeper)) <= length(params.checkpoint_sequence_locs)/2
%                       final = (params.checkpoint_sequence_locs((ceil(pinch_data_count/(100*timeper))-1)*2+1)+params.checkpoint_sequence_locs((ceil(pinch_data_count/(100*timeper))-1)*2+2))/2;
%                     else
%                       final = 0;
%                     end
%                     if ceil(pinch_data_count/(100*timeper)) ==1
%                       initial = 0;
%                     elseif ceil(pinch_data_count/(100*timeper)) <= length(params.checkpoint_sequence_locs)/2
%                       initial = (params.checkpoint_sequence_locs((ceil(pinch_data_count/(100*timeper))-2)*2+1)+params.checkpoint_sequence_locs((ceil(pinch_data_count/(100*timeper))-2)*2+2))/2;
%                     else
%                       initial = 0;
%                     end
%                     
%                     pinch_data_temp = final + (initial-final)*exp(-(rem(pinch_data_count,100*timeper))/100/0.5);


                    current_function = random_function(trialcount,1);
                    current_a = random_function(trialcount,2);
                    current_b = random_function(trialcount,3);
                    current_c = random_function(trialcount,4);
                    current_a2 = random_function(trialcount,5);
                    current_b2 = random_function(trialcount,6);
                    
                    %curr_interleaved_function = interleaved_function(trialcount,1);
                    %curr_interleaved_a = interleaved_function(trialcount,2);
                    %curr_interleaved_b = interleaved_function(trialcount,3);
                    %curr_interleaved_c = interleaved_function(trialcount,4);
                    
                    
                    %now display it
                    %                     if rem(pinch_data_count,10)==0
                    if selected_function ==1 %linear
                        location = [a*pinch_data_temp + c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)];
                    elseif selected_function ==2 %exponential
                        location = [a*b^pinch_data_temp + c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)];
                    elseif selected_function == 3 %logarithm
                        location = [log(a*pinch_data_temp+1)*((a*pinch_data_temp+1)>0)/log(b)+ c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)];
                        if isnan(location)
                            location = 0;
                        end
                    elseif selected_function ==4 %random
                        if current_function ==1 %linear
                            location = [current_a*pinch_data_temp + current_c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)]; 
                        elseif current_function ==2 %exponential
                            location = [current_a*current_b^pinch_data_temp + current_c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)];
                        elseif current_function ==3 %logarithm
                            location = [log(current_a*pinch_data_temp+1)*((current_a*pinch_data_temp+1)>0)/log(current_b)+ current_c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)];
                            if isnan(location)
                                location = 0;
                            end
                        elseif current_function ==4 %sigmoid
                            location = [1./(1+(pinch_data_temp/current_b).^-current_a) + current_c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)]; 
                        else %double sigmoid
                            location = [((1./(1+exp(-current_a*(pinch_data_temp-current_b))))+(1./(1+exp(-current_a2*(pinch_data_temp-current_b2)))))./2 + current_c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)]; 
                        end
                    elseif selected_function == 5 %interleaved
                        if curr_interleaved_function == 1 % logarithmic
                            set(hMainFigure,'Color',get(0, 'defaultuicontrolbackgroundcolor'))
                            % This sets up the initial plot with the bars
%                             for i = 1:size(params.checkpoint_locs,2)
%                                 hCheckpoint(i) = rectangle('Parent',hPlotAxes,'Position',[params.checkpoint_locs(i)-params.checkpoint_thickness(i)/2,params.checkpoint_center_loc(i)-params.checkpoint_height(i)/2,params.checkpoint_thickness(i),params.checkpoint_height(i)],'FaceColor',[params.checkpoint_color_red(i),params.checkpoint_color_green(i),params.checkpoint_color_blue(i)],'LineStyle','none');
%                             end

                            location = [log(curr_interleaved_a*pinch_data_temp+1)*((curr_interleaved_a*pinch_data_temp+1)>0)/log(curr_interleaved_b)+ curr_interleaved_c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)];
                            if isnan(location)
                                location = 0;
                            end
                        elseif curr_interleaved_function == 2 % exponential
                            set(hMainFigure,'Color',[0.2 0.3 0.2])
%                             for i = 1:size(params.checkpoint_locs,2)
%                                 hCheckpoint(i) = rectangle('Parent',hPlotAxes,'Position',[params.checkpoint_locs(i)-params.checkpoint_thickness(i)/2,params.checkpoint_center_loc(i)-params.checkpoint_height(i)/2,params.checkpoint_thickness(i),params.checkpoint_height(i)],'FaceColor',[1 1 1],'EdgeColor',[params.checkpoint_color_red(i),params.checkpoint_color_green(i),params.checkpoint_color_blue(i)],'LineStyle','-','LineWidth',2);
%                             end
                            location = [curr_interleaved_a*curr_interleaved_b^pinch_data_temp + curr_interleaved_c - params.pinchcursor_thickness(1)/2,params.pinchcursor_height_loc-params.pinchcursor_thickness(2)/2,params.pinchcursor_thickness(1),params.pinchcursor_thickness(2)];
                        end
                    end
                    set(hRect,'Position',location);
                    
                    
                    %location is not actually saved, but is stored in
                    %memory...just incase
                    location_data_temp = location(1)+params.pinchcursor_thickness(1)/2;

                    minutes = floor(floor((toc-tp)-tstart)/60);
                    seconds = rem(floor((toc-tp)-tstart),60);
%                     set(hTimerbox,'String',strcat(int2str(minutes),':',int2str(seconds)));
                    set(hTimerbox,'String',(0));
                    %in case needed, use below to display voltage values in
                    %timer box instead of the time
                    %set(hTimerbox,'String',strcat(num2str(pinch_data_temp)));
                
                    
                    %mouseover condition
                    mouseover_checkpoint = find ((location_data_temp<(params.checkpoint_locs+params.checkpoint_thickness/2)).*(location_data_temp>(params.checkpoint_locs - params.checkpoint_thickness/2)));
                    if mouseover == 0 && ~isempty(mouseover_checkpoint)
                    set(hCheckpoint(mouseover_checkpoint),'FaceColor',[params.mouseover_checkpoint_color_red(mouseover_checkpoint),params.mouseover_checkpoint_color_green(mouseover_checkpoint),params.mouseover_checkpoint_color_blue(mouseover_checkpoint)]);
                    mouseover = mouseover_checkpoint;
                    elseif mouseover ~= 0 && isempty(mouseover_checkpoint)
                        set(hCheckpoint(mouseover),'FaceColor',[params.checkpoint_color_red(mouseover),params.checkpoint_color_green(mouseover),params.checkpoint_color_blue(mouseover)]);
                        mouseover = 0;
                    end
                        

                    %stopcondition
                    if pinch_data_count >(params.sample_rate*params.trialwait)
                        break;
                    else %follow the bars algorithm
                        if (location_data_temp<params.checkpoint_sequence_locs((cursor_loc-1)*2+2))&&...
                                (location_data_temp>params.checkpoint_sequence_locs((cursor_loc-1)*2+1))&&...
                                stays_in==0

                            stays_in = 1;
                            time_start = (toc-tp);
                        elseif  (location_data_temp<params.checkpoint_sequence_locs((cursor_loc-1)*2+2))&&...
                                (location_data_temp>params.checkpoint_sequence_locs((cursor_loc-1)*2+1))&&...
                                (stays_in==1)&&((toc-tp)-time_start)>=params.checkpoint_holdtime(cursor_loc)
                            stays_in = 0; %move to the next checkpoint
                            cursor_loc = cursor_loc+1;
                            if cursor_loc > length(params.checkpoint_sequence_locs)/2
                                set(hTrialDone,'FaceColor','r');
                               
                                break %done!!!
                            end
                        elseif (location_data_temp<params.checkpoint_sequence_locs((cursor_loc-1)*2+2))&&...
                                (location_data_temp>params.checkpoint_sequence_locs((cursor_loc-1)*2+1))&&...
                                (stays_in==1)
                            %keep going!
                        else
                            stays_in = 0; %gotta finish it!!
                        end
                    end
                    
                    drawnow                    
                    
                    
                    %increment count
                    pinch_data_count = pinch_data_count + 1;
                    %display time
                    e(pinch_data_count) = (toc-tp);
                    %record voltage
                    pinch_data(pinch_data_count,1) = pinch_data_temp;
                    location_data(pinch_data_count) = location_data_temp;
                    
                    if (pinch_data_count>1) && ((e(pinch_data_count) - e(pinch_data_count-1)) > 1/params.sample_rate)
                        lag_count = lag_count + 1;
                    end
                    

%                     %calculate how much time needs to be passed
%                     next_time = floor((toc-tp)/params.sample_rate)+ 1/params.sample_rate; %initialization time
%                     %calculate how much time has passed by
%                     while ((toc-tp)<next_time)
%                         %wait
%                     end

                end
                %take down end time
                tend = (toc-tp);
                %end of trial
                TrialStartButton_status = 0;
                set(hTrialStartButton,'String','Pause');

                
                
                %save raw data
                fprintf(fid_raw,'%-3.0f\n',trialcount); %note how x and y are placed
                %print time (number of seconds since start of first trial)
                fprintf(fid_raw,'%-10.10f\n',tstart);
                fprintf(fid_raw,'%-10.10f\n',tend);
                fprintf(fid_raw,'%-10.10f\n',lag_count);
                
                size(pinch_data)
                size(e)
                pinch_data_count
                for i = 1:pinch_data_count
                      fprintf(fid_raw,'%10.10f\t%10.10f\t%10.10f\t%10.10f\t%1.3f\n',[pinch_data(i).*params.internal_calibration_gain;pinch_data(i).*params.internal_calibration_gain*-129+158;e(i);location_data(i);TMStrigger_status])
                end  
                
                
                %process raw information in order to extract perfect 100 Hz
                %signal that can be used for analysis
                e = e-e(1);
                
                %remove all duplicate time points (shouldnt occur)
                for i = 1:pinch_data_count
                    if sum(e==e(i))>1
                        temp_loc = find(e(i)==e);
                        temp_loc(1) = [];
                        e(temp_loc) = [];
                        pinch_data(temp_loc) = [];
                        location_data(temp_loc) = [];
                    end
                end
                
                %process interpolation algorithm
                t = 0;
                i = 1;
                while(sum(t<=e)>0)
                    loc_pinch_data_2 = find(t<=e,1,'first');
                    loc_pinch_data_1 = find(t>=e,1,'last');
                    if ( loc_pinch_data_1~= loc_pinch_data_2)
                        pinch_data_clean(i) = (t-e(loc_pinch_data_1))/(e(loc_pinch_data_2)-e(loc_pinch_data_1))*...
                            (pinch_data(loc_pinch_data_2)-pinch_data(loc_pinch_data_1)) + pinch_data(loc_pinch_data_1);
                        location_data_clean(i) = (t-e(loc_pinch_data_1))/(e(loc_pinch_data_2)-e(loc_pinch_data_1))*...
                            (location_data(loc_pinch_data_2)-location_data(loc_pinch_data_1)) + location_data(loc_pinch_data_1);
                        e_clean(i) = (t-e(loc_pinch_data_1))/(e(loc_pinch_data_2)-e(loc_pinch_data_1))*...
                            (e(loc_pinch_data_2)-e(loc_pinch_data_1)) + e(loc_pinch_data_1);
                    else
                        pinch_data_clean(i) = pinch_data(loc_pinch_data_1);
                        location_data_clean(i) = location_data(loc_pinch_data_1);
                        e_clean(i) = e(loc_pinch_data_1);
                    end                    
                    t = t+ 1/params.sample_rate;
                    i = i+1;
                end
               %happy = [pinch_data_clean(i).*params.internal_calibration_gain;pinch_data_clean(i).*params.internal_calibration_gain*-129+158;TMStrigger_status]

                                %turn lag box on
%                 set(hLagbox,'Visible','on','String',num2str(lag_count/pinch_data_count * 100));
%                 set(hLagText,'Visible','on');
                
                

                %save data
                %print trial in file
                fprintf(fid,'%-3.0f\n',trialcount); %note how x and y are placed

                %print time (number of seconds since start of first trial)
                fprintf(fid,'%-10.0f\n',floor(tend*1000));

                %save data calibrated to the internal gain set the way
                %kopal set it in her labview program
                for i = 1:length(e_clean)
                      fprintf(fid,'%1.3f\t%1.3f\t%1.3f\n',[pinch_data_clean(i).*params.internal_calibration_gain;pinch_data_clean(i).*params.internal_calibration_gain*-129+158;TMStrigger_status]);
%                       fprintf(fid_raw,'%1.3f\t%1.3f\t%1.3f\n',[location(2);pinch_data_clean(i).*params.internal_calibration_gain*-129+158;pinch_data_clean(i)])

                end
                
                %reset variables
                e = [];
                pinch_data = [];
                location_data = [];
                e_clean = [];
                pinch_data_clean = [];
                location_data_clean = [];

                if (StartButton_status == 0)
                    fclose('all')
                    break;
                end
                
                if params.manual_next_trial == 1
                    %set(hTrialStartButton,'Visible','on')
                    %beep;
                    
                    
                    %check for program stop or next trial stop
                    %uiwait;
                    %pause(1);
                    sound(beep1);
                    disp('2')
                    pause(1);
                    sound(beep1);
                    set(hTrialDone,'FaceColor','y');
                    disp('1')
                    pause(1);
                    sound(beep2);
                    set(hTrialDone,'FaceColor','g');
                    disp('GO!')
                    
                    %set(hTrialStartButton,'Visible','off')
                    trialcount;
                    if (StartButton_status == 0)
                        fclose('all')
                        break;
                    end
                    tdelay = (toc-tp);
                    while (tdelay-tend)<params.delay_next_trial
                        tdelay = (toc-tp);
                    end
                else
                    tdelay = (toc-tp);
                    while (tdelay-tend)<params.delay_next_trial
                        tdelay = (toc-tp);
                    end
                end

                if (StartButton_status == 0)
                    fclose('all')
                    break;
                end
                
            end
            fclose('all')
            
            tend = (toc-tp);
            minutes = floor(floor((toc-tp)-tstart)/60);
            seconds = rem(floor((toc-tp)-tstart),60);
            set(hTimerbox,'String',strcat(int2str(minutes),':',int2str(seconds)));

            
            
            StartButton_status = 0;
            set(hStartButton,'String','START');
            uiresume;
            hbutt = rectangle('Parent',hReadyAxes,'Position',[-0.5,-0.5,1,1],'FaceColor','r','Curvature',1);
            %set(hTrialStartButton,'Visible','off')



            %make the boxes invisible
            set(hFunctionPopupmenu,'Visible','on');
            set(hGain1box,'Visible','on');
            set(hGain2box,'Visible','on');
            set(hGain3box,'Visible','on');
            set(hTMScheck,'Visible','on');
            set(hMRIcheck,'Visible','on');
            set(hFunctionPopupmenuText,'Visible','on');
            set(hGain1boxText,'Visible','on');
            set(hGain2boxText,'Visible','on');
            set(hGain3boxText,'Visible','on');

            set(hLagbox,'Visible','on');
            set(hLagText,'Visible','on');

            set(hTrialbox,'String',num2str(ntrials));
            set(hFunctionPopupmenu,'value',selected_function);
            set(hGain1box,'String',num2str(a));
            set(hGain2box,'String',num2str(b));
            set(hGain3box,'String',num2str(c));
            set(hTMScheck,'String',tms_trigger);
            set(hMRIcheck,'String',mri_trigger);

        else
            StartButton_status = 0;
            set(hStartButton,'String','START');
            uiresume;
            hbutt = rectangle('Parent',hReadyAxes,'Position',[-0.5,-0.5,1,1],'FaceColor','r','Curvature',1);
            %set(hTrialStartButton,'Visible','off')

            %make the boxes invisible
            set(hFunctionPopupmenu,'Visible','on');
            set(hGain1box,'Visible','on');
            set(hGain2box,'Visible','on');
            set(hGain3box,'Visible','on');
            set(hTMScheck,'Visible','on');
            set(hMRIcheck,'Visible','on');
            set(hFunctionPopupmenuText,'Visible','on');
            set(hGain1boxText,'Visible','on');
            set(hGain2boxText,'Visible','on');
            set(hGain3boxText,'Visible','on');

            set(hTrialbox,'String',num2str(ntrials));
            set(hFunctionPopupmenu,'value',selected_function);
            set(hGain1box,'String',num2str(a));
            set(hGain2box,'String',num2str(b));
            set(hGain3box,'String',num2str(c));
            set(hTMScheck,'String',tms_trigger);
            set(hMRIcheck,'String',mri_trigger);
        end

    end

    function nextTrial
        TrialStartButton_status = 1;
        uiresume;
        set(hTrialStartButton,'Visible','off')
    end
%----------------------------------------------------------------------
    function hTrialStartButtonCallback(hObject, eventdata)
        if StartButton_status == 1
            toc
            t1=tic
            if TrialStartButton_status == 1
                set(hTrialStartButton,'String','Resume');
                disp(TrialStartButton_status)
                TrialStartButton_status = 0
                uiwait;
            else
                set(hTrialStartButton,'String','Pause');
                disp(TrialStartButton_status)
                TrialStartButton_status = 1
                uiresume;
            end
            tp=tp+toc(t1)
        end
        %uiresume;
        %set(hTrialStartButton,'Visible','off')
    end

%----------------------------------------------------------------------
    function hFileboxCallback(hObject, eventdata)
        % Callback function run when the update button is pressed
        %         localUpdatePlot();
    end

%----------------------------------------------------------------------
    function hMRIcheckCallback(hObject, eventdata)
        % Callback function run when the update button is pressed
        %         localUpdatePlot();
        status = get(hMRIcheck,'value');
        MRItrigger_status = status;
    end

%----------------------------------------------------------------------
    function hTMScheckCallback(hObject, eventdata)
        % Callback function run when the update button is pressed
        %         localUpdatePlot();
        status = get(hTMScheck,'value');
        TMStrigger_status = status;
    end

%----------------------------------------------------------------------
    function hOpenMenuitemCallback(hObject, eventdata)
        % Callback function run when the Open menu item is selected
        file = uigetfile('*.fig');
        if ~isequal(file, 0)
            open(file);
        end
    end

%----------------------------------------------------------------------
    function hPrintMenuitemCallback(hObject, eventdata)
        % Callback function run when the Print menu item is selected
        printdlg(hMainFigure);
    end

%----------------------------------------------------------------------
    function hCloseMenuitemCallback(hObject, eventdata)
        % Callback function run when the Close menu item is selected
        selection = questdlg(['Close ' get(hMainFigure,'Name') '?'],...
            ['Close ' get(hMainFigure,'Name') '...'],...
            'Yes','No','Yes');
        if strcmp(selection,'No')
            return;
        end

        delete(hMainFigure);
    end


end % end of axesMenuToolbar
