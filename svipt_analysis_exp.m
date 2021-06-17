
tic
clear all;
fclose all;
close all;

%% Parameter to process Force %%
fs_Pos = 100; %Sampling rate of JH (since MEP of TMS can be chopped with fs of 1000)
T_Pos = 1 / fs_Pos;
fn_Pos = fs_Pos / 2;
fre1_Pos = 10;
fre2_Pos = 20;
fre3_Pos = 6;
fre4_Pos = 12;
f1n_Pos = fre1_Pos / fn_Pos;
f2n_Pos = fre2_Pos / fn_Pos;
f3n_Pos = fre3_Pos / fn_Pos;
f4n_Pos = fre4_Pos / fn_Pos;

[Pos_b, Pos_a]=butter(4,f2n_Pos,'low');
[Pos_d, Pos_c]=butter(4,f4n_Pos,'low');

%% Importing file with information for processing
%StartDirectory = cd('/Volumes/RONAN_USB/Celnik_Postdoc/Experiments/NML_EI/Data/Behavior/x_New_SVIPT_analysis_x');
StartDirectory = cd('/Volumes/RONAN_USB/Celnik_Postdoc/Experiments/NML_EI/Data/Behavior/SVIPT_analysis');
Subjects_Full = importdata('Test.txt');
%Subject =  Subjects_Full(1:2:end);
Subject =  Subjects_Full(1);

% Day_num = ['1'; '2'];
% Block_num = ['1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9';];

% Day_num = '1';
% Block_num = '1';

% EMG_Folder = 'TMS';
Force_Folder = 'raw_int';

%% SVIPT INFO
gate2_low = 0.169-0.0155;
gate2_high = 0.264+0.0155;
gate2_mid = ((gate2_high-gate2_low)/2)+gate2_low;

gate4_low = 0.359-0.0155;
gate4_high = 0.435+0.0155;
gate4_mid = ((gate4_high-gate4_low)/2)+gate4_low;

gate1_low = 0.551-0.0155;
gate1_high = 0.627+0.0155;
gate1_mid = ((gate1_high-gate1_low)/2)+gate1_low;

gate5_low = 0.715-0.03;
gate5_high = 0.715+0.03;

gate3_low = 0.831-0.0155;
gate3_high = 0.929+0.0155;
gate3_mid = ((gate3_high-gate3_low)/2)+gate3_low;

%% ANALYSIS STARTS %%
%Save_dir = sprintf ('//Users//agostinacasamentomoran//Documents//Documents//JH_Celnik_Lab//Projects//1.MotorLearning_Fatigue_MMActivity//EMG_ANALYSIS_OUTPUT');
%Save_dir = sprintf ('/Volumes/RONAN_USB/Celnik_Postdoc/Experiments/NML_EI/Data/Behavior/x_New_SVIPT_analysis_x/Output'); % Save location for output file
Save_dir = sprintf ('/Volumes/RONAN_USB/Celnik_Postdoc/Experiments/NML_EI/Data/Behavior/SVIPT_analysis/Output'); % Save location for output file

for s = 1:length(Subject) % Getting all participants
    Folder_Name = [Subject{s}]; %folder name
    Participant_OutputFolder_Name = sprintf('%s//%s_Analysis',Save_dir, Folder_Name);
    mkdir (Participant_OutputFolder_Name);      % creates a new folder for participant folder in Analysis Output folder
    
    % loading FORCE data from "block" folder
    %DirectoryPathName2 = sprintf ('//Users//agostinacasamentomoran//Documents//Documents//JH_Celnik_Lab//Projects//1.MotorLearning_Fatigue_MMActivity//RAW_DATA//%s//d%s//%s',Folder_Name,Day_num(d),Force_Folder); %Name of the directory
    %DirectoryPathName2 = sprintf ('/Volumes/RONAN_USB/Celnik_Postdoc/Experiments/NML_EI/Data/Behavior/x_New_SVIPT_analysis_x/Data/%s/%s',Folder_Name,Force_Folder);
    DirectoryPathName2 = sprintf ('/Volumes/RONAN_USB/Celnik_Postdoc/Experiments/NML_EI/Data/Behavior/SVIPT_analysis/Data/%s/%s',Folder_Name,Force_Folder);
    WorkingDirectory2 = cd (DirectoryPathName2);
    
    [FileList, ~] = glob(['*raw']);
    ArrangedFileList_Force = sort_nat(FileList);
    
    
    %% Initializing the stats matrices
    all_day_avg = []; %Initializing the matrix that will contain the avg of each block
    all_day_sd = []; %Initializing the matrix that will contain the sd of each block
    all_day_cv = []; %Initializing the matrix that will contain the cv of each block
    all_day_sem = []; %Initializing the matrix that will contain the sem of each block
    
    %% Creating Output file for stats
    cd(Participant_OutputFolder_Name)
    outfile_AVG = sprintf('%s_blocks_SVIPTAnalysis_AVG.txt',Folder_Name);
    fid_AVG = fopen(outfile_AVG,'w');
    fprintf(fid_AVG,'Block_Num\t Original_Pk_Ct\t Bin_error_all\t Bin_error_percent\t abs_Dist_error_center\t Dist_error_center\t abs_Dist_error_gates\t Dist_error_gates\t trial_duration\t Bin_error_pk1\t Bin_error_pk2\t Bin_error_pk3\t Bin_error_pk4\t');
    fprintf(fid_AVG,'abs_Dist_error_ce_pk1\t abs_Dist_error_ce_pk2\t abs_Dist_error_ce_pk3\t abs_Dist_error_ce_pk4\t Dist_error_ce_pk1\t Dist_error_ce_pk2\t Dist_error_ce_pk3\t Dist_error_ce_pk4\t'); 
    fprintf(fid_AVG,'abs_Dist_error_ga_pk1\t abs_Dist_error_ga_pk2\t abs_Dist_error_ga_pk3\t abs_Dist_error_ga_pk4\t Dist_error_ga_pk1\t Dist_error_ga_pk2\t Dist_error_ga_pk3\t Dist_error_ga_pk4\t ttpeak1\t ttpeak2\t ttpeak3\t ttpeak4\t ttpeak5\t');
    fprintf(fid_AVG,'ttrelax1\t ttrelax2\t ttrelax3\t ttrelax4\t ttrelax5\t duration1\t duration2\t duration3\t duration4\t duration5\t pk_force1\t pk_force2\t pk_force3\t pk_force4\t pk_force5\n');
% 
% outfile_SD = sprintf('%s_d%s_blocks_SVIPTAnalysis_SD.txt',Folder_Name,Day_num(d));
% fid_SD = fopen(outfile_SD,'w');
% fprintf(fid_SD,'Block_Num\t Original_Pk_Ct\t Bin_error_all\t Bin_error_percent\t abs_Dist_error_center\t Dist_error_center\t abs_Dist_error_gates\t Dist_error_gates\t trial_duration\t Bin_error_pk1\t Bin_error_pk2\t Bin_error_pk3\t Bin_error_pk4\t');
% fprintf(fid_SD,'abs_Dist_error_ce_pk1\t abs_Dist_error_ce_pk2\t abs_Dist_error_ce_pk3\t abs_Dist_error_ce_pk4\t Dist_error_ce_pk1\t Dist_error_ce_pk2\t Dist_error_ce_pk3\t Dist_error_ce_pk4\t');
% fprintf(fid_SD,'abs_Dist_error_ga_pk1\t abs_Dist_error_ga_pk2\t abs_Dist_error_ga_pk3\t abs_Dist_error_ga_pk4\t Dist_error_ga_pk1\t Dist_error_ga_pk2\t Dist_error_ga_pk3\t Dist_error_ga_pk4\t ttpeak1\t ttpeak2\t ttpeak3\t ttpeak4\t ttpeak5\t');
% fprintf(fid_SD,'ttrelax1\t ttrelax2\t ttrelax3\t ttrelax4\t ttrelax5\t duration1\t duration2\t duration3\t duration4\t duration5\t pk_force1\t pk_force2\t pk_force3\t pk_force4\t pk_force5\n');
% 
% outfile_CV = sprintf('%s_d%s_blocks_SVIPTAnalysis_CV.txt',Folder_Name,Day_num(d));
% fid_CV = fopen(outfile_CV,'w');
% fprintf(fid_CV,'Block_Num\t Original_Pk_Ct\t Bin_error_all\t Bin_error_percent\t abs_Dist_error_center\t Dist_error_center\t abs_Dist_error_gates\t Dist_error_gates\t trial_duration\t Bin_error_pk1\t Bin_error_pk2\t Bin_error_pk3\t Bin_error_pk4\t');
% fprintf(fid_CV,'abs_Dist_error_ce_pk1\t abs_Dist_error_ce_pk2\t abs_Dist_error_ce_pk3\t abs_Dist_error_ce_pk4\t Dist_error_ce_pk1\t Dist_error_ce_pk2\t Dist_error_ce_pk3\t Dist_error_ce_pk4\t');
% fprintf(fid_CV,'abs_Dist_error_ga_pk1\t abs_Dist_error_ga_pk2\t abs_Dist_error_ga_pk3\t abs_Dist_error_ga_pk4\t Dist_error_ga_pk1\t Dist_error_ga_pk2\t Dist_error_ga_pk3\t Dist_error_ga_pk4\t ttpeak1\t ttpeak2\t ttpeak3\t ttpeak4\t ttpeak5\t');
% fprintf(fid_CV,'ttrelax1\t ttrelax2\t ttrelax3\t ttrelax4\t ttrelax5\t duration1\t duration2\t duration3\t duration4\t duration5\t pk_force1\t pk_force2\t pk_force3\t pk_force4\t pk_force5\n');

% outfile_SEM = sprintf('%s_d%s_blocks_SVIPTAnalysis_SEM.txt',Folder_Name,Day_num(d));
% fid_SEM = fopen(outfile_SEM,'w');
% fprintf(fid_SEM,'Block_Num\t Original_Pk_Ct\t Bin_error_all\t Bin_error_percent\t abs_Dist_error_center\t Dist_error_center\t abs_Dist_error_gates\t Dist_error_gates\t trial_duration\t Bin_error_pk1\t Bin_error_pk2\t Bin_error_pk3\t Bin_error_pk4\t');
% fprintf(fid_SEM,'abs_Dist_error_ce_pk1\t abs_Dist_error_ce_pk2\t abs_Dist_error_ce_pk3\t abs_Dist_error_ce_pk4\t Dist_error_ce_pk1\t Dist_error_ce_pk2\t Dist_error_ce_pk3\t Dist_error_ce_pk4\t');
% fprintf(fid_SEM,'abs_Dist_error_ga_pk1\t abs_Dist_error_ga_pk2\t abs_Dist_error_ga_pk3\t abs_Dist_error_ga_pk4\t Dist_error_ga_pk1\t Dist_error_ga_pk2\t Dist_error_ga_pk3\t Dist_error_ga_pk4\t ttpeak1\t ttpeak2\t ttpeak3\t ttpeak4\t ttpeak5\t');
% fprintf(fid_SEM,'ttrelax1\t ttrelax2\t ttrelax3\t ttrelax4\t ttrelax5\t duration1\t duration2\t duration3\t duration4\t duration5\t pk_force1\t pk_force2\t pk_force3\t pk_force4\t pk_force5\n');

    for b = 1:length(ArrangedFileList_Force) % Looping through all identified *raw files
        file = ArrangedFileList_Force{b};
        fileData = dlmread(file,'');
        
        Trial_info = fileData(:,1); % Importing column with trial #,tstart, tend, logcount
        Time_pts = fileData(:,3); % Importing time
        Force_raw = fileData(:,4); % Importing force
        
        Force = filtfilt(Pos_b, Pos_a, Force_raw);
        
        all_block_data = []; % Initializing matrix that will contrain the data for all trials of each block
        fig1_num = b;
        fig2_num = fig1_num + 1;
        
        fig1_name = num2str(fig1_num);
        fig2_name = num2str(fig2_num);
        
        %% Creating Output file
        cd(Participant_OutputFolder_Name)
        outfile1 = sprintf('%s_b%s_Analysis.txt',Folder_Name,fig1_name);
        fid1 = fopen(outfile1,'w');
        fprintf(fid1,'Trial_Num\t Original_Pk_Ct\t Bin_error_all\t Bin_error_percent\t abs_Dist_error_center\t Dist_error_center\t abs_Dist_error_gates\t Dist_error_gates\t trial_duration\t Bin_error_pk1\t Bin_error_pk2\t Bin_error_pk3\t Bin_error_pk4\t');
        fprintf(fid1,'abs_Dist_error_ce_pk1\t abs_Dist_error_ce_pk2\t abs_Dist_error_ce_pk3\t abs_Dist_error_ce_pk4\t Dist_error_ce_pk1\t Dist_error_ce_pk2\t Dist_error_ce_pk3\t Dist_error_ce_pk4\t');
        fprintf(fid1,'abs_Dist_error_ga_pk1\t abs_Dist_error_ga_pk2\t abs_Dist_error_ga_pk3\t abs_Dist_error_ga_pk4\t Dist_error_ga_pk1\t Dist_error_ga_pk2\t Dist_error_ga_pk3\t Dist_error_ga_pk4\t ttpeak1\t ttpeak2\t ttpeak3\t ttpeak4\t ttpeak5\t');
        fprintf(fid1,'ttrelax1\t ttrelax2\t ttrelax3\t ttrelax4\t ttrelax5\t duration1\t duration2\t duration3\t duration4\t duration5\t pk_force1\t pk_force2\t pk_force3\t pk_force4\t pk_force5\n');

        %% Quantifying some properties of each practice block
        Real_duration_sec = Time_pts(end); % getting block duration
        Real_duration_5000Hz = Real_duration_sec * 5000; % converting block duration into 5000 Hz
        
        [pks_y,pks_x] = findpeaks(Force,'MinPeakProminence',0.15); % finding pk forces with drops higher than 0.1 in b1
        
        % Getting each trial start & end time
        zeroes = Time_pts == 0; % ID rows that actually contain trial #,tstart, tend, logcount
        
        % Getting rid of lagcount
        [indeces, values] = find(zeroes == 1);
        lagcount_indeces = indeces(4:4:end);
        zeroes(lagcount_indeces) = 0;
        
        Trial_info_filt = Trial_info(zeroes); % Getting only trial num, start and end of trial (not lag count)
        numb_trials = Trial_info_filt(1:3:end-2); % Getting number of trials performed
                
        tstart_sec = Trial_info_filt(2:3:end-1);% Getting only start time in sec
        tstart_100Hz = tstart_sec * 100; % Converting start time to 100Hz
        
        tend_sec = Trial_info_filt(3:3:end);% Getting only end time in sec
        tend_100Hz = tend_sec * 100; % Converting end time to 100Hz
        
        tduration_sec =  tend_sec - tstart_sec; % Getting trial duration in sec
        tduration_100Hz =  tduration_sec * 100; % Converting trial duration to 100Hz
        
        %find trial break points in force data
        Trial_info_nolag = Trial_info .* (double(zeroes));
        breakpoints_force = [];
        for trl = 1: length(numb_trials)
            bkpt = find(Trial_info_nolag == trl, 1, 'first');
            breakpoints_force = [breakpoints_force, bkpt];
        end
        breakpoints_force =  [breakpoints_force,length(Force)];
        
        
        %% Getting & plotting each force trial data
        % Initializing all_trial matrix
        all_force_trls = zeros([round(max(tduration_100Hz)) length(tduration_100Hz)]); % Col = # of ID trls; Row = length of longest ID trl
        
        trial_Count = [1:(length(breakpoints_force)-1)];
        
        for ft = 1:(length(breakpoints_force)-1) % Trial loop starts
            trl_num = num2str(ft);
            force_trial_raw = Force_raw(breakpoints_force(ft):breakpoints_force(ft+1));
            force_trial = Force(breakpoints_force(ft):breakpoints_force(ft+1));
            
            force_trial_time = Time_pts(breakpoints_force(ft):breakpoints_force(ft+1));
            
            [trl_pks_y,trl_pks_x] = findpeaks(force_trial,'MinPeakProminence',0.15); % finding pk forces with drops higher than 0.1 in b1
            original_pk_count = length(trl_pks_y);
            original_pks = trl_pks_y;
            
            exclude_trial = 0;
            
            if (length(trl_pks_y) < 5) || (length(trl_pks_y) > 5)
                
                fig_title = sprintf('Plotting Block %s Trial %s', fig1_name, trl_num);
                figure(200)
                plot(force_trial)
                title(fig_title)
                hold on
                yline(gate1_low,'-b');
                yline(gate1_high,'-b','Gate 1', 'LabelVerticalAlignment', 'middle');
                
                yline(gate2_low,'-r');
                yline(gate2_high,'-r','Gate 2', 'LabelVerticalAlignment', 'middle');
                
                yline(gate3_low,'-g');
                yline(gate3_high,'-g','Gate 3', 'LabelVerticalAlignment', 'middle');
                
                yline(gate4_low,'-','color', '#D95319');
                yline(gate4_high,'-','Gate 4', 'LabelVerticalAlignment', 'middle','color', '#D95319');
                
                yline(gate5_low,'-m');
                yline(gate5_high,'-m','Gate 5', 'LabelVerticalAlignment', 'middle');
                hold off
                
                abort = menu('Do you want to exclude this trial? Press Y or N', 'Y', 'N');
                if abort == 1
                    close (figure(200))
                    exclude_trial = 1;
                elseif abort == 2
                    close (figure(200))
                    parameter = 1;
                    while parameter == 1
                        fig_title = sprintf('Plotting Block %s Trial %s', fig1_name, trl_num);
                        
                        figure(100)
                        plot(force_trial)
                        title(fig_title)
                        hold on
                        yline(gate1_low,'-b');
                        yline(gate1_high,'-b','Gate 1', 'LabelVerticalAlignment', 'middle');
                        
                        yline(gate2_low,'-r');
                        yline(gate2_high,'-r','Gate 2', 'LabelVerticalAlignment', 'middle');
                        
                        yline(gate3_low,'-g');
                        yline(gate3_high,'-g','Gate 3', 'LabelVerticalAlignment', 'middle');
                        
                        yline(gate4_low,'-','color', '#D95319');
                        yline(gate4_high,'-','Gate 4', 'LabelVerticalAlignment', 'middle','color', '#D95319');
                        
                        yline(gate5_low,'-m');
                        yline(gate5_high,'-m','Gate 5', 'LabelVerticalAlignment', 'middle');
                        
                        [trl_pks_x,trl_pks_y] = ginput(5); % ask for user input to get peaks
                        plot(force_trial,'k')
                        title(fig_title)
                        plot(trl_pks_x,trl_pks_y,':k*')
                        plot(trl_pks_x(1),trl_pks_y(1),'b*')
                        plot(trl_pks_x(2),trl_pks_y(2),'r*')
                        plot(trl_pks_x(3),trl_pks_y(3),'g*')
                        plot(trl_pks_x(4),trl_pks_y(4),'*','color', '#D95319')
                        plot(trl_pks_x(5),trl_pks_y(5),'m*')
                        hold off
                        
                        
                        choice = menu('Is your selection correct? Press Y or N', 'Y', 'N');
                        if choice ~= 1
                            parameter = 1;
                        else
                            parameter = 0;
                        end
                        close (figure(100))
                    end
                end
            end
            
            selected_peaks = trl_pks_y;
            
            % Analyze only the trials that we want to include, make all outcome variable zeros for the trials that we want to exlude
            if exclude_trial ~= 1
                %% Quantified Binary Error %%
                Bin_error = [0, 0, 0, 0]; %Initializing bin error assuming no error for all 4 peaks
                all_low = [gate1_low, gate2_low, gate3_low, gate4_low];
                all_high = [gate1_high, gate2_high, gate3_high, gate4_high];
                for bn = 1:length(Bin_error)
                    if trl_pks_y(bn) < all_low(bn) || trl_pks_y(bn) > all_high(bn) %Deciding if peak outside gate
                        Bin_error(bn) = 1; %If outside bin error for each peak changes to 1
                    end
                end
                Bin_error_count = sum(Bin_error); %Total number of bin errors within trial
                if Bin_error_count > 0
                    Bin_error_all = 1; %At least 1 bin error within trial
                    Bin_error_percent = (Bin_error_count/4)*100; %Percentage of peaks with bin error
                else
                    Bin_error_all = 0; %No bin error within trial
                    Bin_error_percent = 0; %No percentage of peaks with bin error
                end
                
                %% Quantified Distance Error relative to center %%
                abs_Dist_error_center = [0, 0, 0, 0]; %Initializing abs dist error assuming no error for all 4 peaks
                Dist_error_center = [0, 0, 0, 0]; %Initializing dist error assuming no error for all 4 peaks
                all_mid = [gate1_mid, gate2_mid, gate3_mid, gate4_mid];
                for ce = 1:length(abs_Dist_error_center)
                    if trl_pks_y(ce) < all_mid(ce) || trl_pks_y(ce) > all_mid(ce)
                        abs_Dist_error_center(ce) = abs(trl_pks_y(ce)-all_mid(ce));
                        Dist_error_center(ce) = trl_pks_y(ce)-all_mid(ce);
                    end
                end
                total_abs_Dist_error_center = sum(abs_Dist_error_center);
                total_Dist_error_center = sum(Dist_error_center);
                
                %% Quantified Distance Error relative to gates %%
                abs_Dist_error_gates = [0, 0, 0, 0]; %Initializing abs dist error assuming no error for all 4 peaks
                Dist_error_gates = [0, 0, 0, 0]; %Initializing dist error assuming no error for all 4 peaks
                for ce = 1:length(abs_Dist_error_gates)
                    if trl_pks_y(ce) < all_low(ce) 
                        abs_Dist_error_gates(ce) = abs(trl_pks_y(ce)-all_low(ce));
                        Dist_error_gates(ce) = trl_pks_y(ce)-all_low(ce);
                    elseif trl_pks_y(ce) > all_high(ce) %Deciding if peak outside gates
                        abs_Dist_error_gates(ce) = abs(trl_pks_y(ce)-all_high(ce));
                        Dist_error_gates(ce) = trl_pks_y(ce)-all_high(ce);
                    end
                end
                total_abs_Dist_error_gates = sum(abs_Dist_error_gates);
                total_Dist_error_gates = sum(Dist_error_gates);
                
                %% Quantified Temporal Variables %%
                onsets = [0, 0, 0, 0, 0];
                offsets = [0, 0, 0, 0, 0];
                ttpeak = [0, 0, 0, 0, 0];
                ttrelax = [0, 0, 0, 0, 0];
                duration = [0, 0, 0, 0, 0];
                
                lowest_pk = min(selected_peaks);
                onset_threshold = lowest_pk * 0.05;
                
                for pk = 1:length(onsets)
                    temp_pk_array_half1 = force_trial_raw(1: trl_pks_x(pk));
                    flip_temp_pk_array = flipud(temp_pk_array_half1);
                    temp_onset = find(flip_temp_pk_array <= onset_threshold, 1, 'first');
                    act_onset = length(temp_pk_array_half1) - temp_onset;
                    act_onset_time = force_trial_time(act_onset);
                    if pk == 1
                        if act_onset_time == 0
                        act_onset_time = force_trial_time(5);
                        end
                    end
                    onsets(pk) = act_onset_time;
                    sample_peak = round(trl_pks_x(pk));
                    time_peak = force_trial_time(sample_peak);
                    act_ttpeak = time_peak - act_onset_time;
                    ttpeak(pk) = act_ttpeak;
                    
                    
                    
                    temp_pk_array_half2 = force_trial_raw(trl_pks_x(pk):end);
                    temp_offset = find(temp_pk_array_half2 <= onset_threshold, 1, 'first');
                    act_offset = length(temp_pk_array_half1) + temp_offset;
                    act_offset_time = force_trial_time(act_offset);
                    offsets(pk) = act_offset_time;
                    act_ttrelax = act_offset_time - time_peak;
                    ttrelax(pk) = act_ttrelax;
                    
                    duration(pk) = act_ttpeak+act_ttrelax;
                    
                end
                
                final_sample_peak = round(trl_pks_x(5));
                final_time_peak = force_trial_time(final_sample_peak);
                trial_duration = final_time_peak-onsets(1);
                
                all_force_trls(1:length(force_trial),ft) = force_trial;
                
                %% Saving all trial variables
                all_trial_data = [ft, original_pk_count, Bin_error_all, Bin_error_percent, total_abs_Dist_error_center, total_Dist_error_center, total_abs_Dist_error_gates, total_Dist_error_gates, trial_duration, ...
                    Bin_error, abs_Dist_error_center, Dist_error_center, abs_Dist_error_gates, Dist_error_gates, ttpeak, ttrelax, duration, trl_pks_y'];
                fprintf(fid1, '%.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n', all_trial_data);
                
                all_block_data = [all_block_data; all_trial_data];
            end
        end
        
        avg_block_data = [];
        sd_block_data = [];
        cv_block_data = [];
        sem_block_data = [];
        
        for var = 2:length(all_block_data) % Quantifying block stats
            trial_count = max(all_block_data(:,1));
            avg_var = mean(all_block_data(:,var));
            sd_var = std(all_block_data(:,var));
            cv_var = sd_var/avg_var*100;
            sem_var = sd_var/sqrt(trial_count);
            
            avg_block_data = [avg_block_data, avg_var];
            sd_block_data = [sd_block_data, sd_var];
            cv_block_data = [cv_block_data, cv_var];
            sem_block_data = [sem_block_data, sem_var];
        end
        
        all_day_avg = [all_day_avg; avg_block_data]; % Adding each block avg in each row
        all_day_sd = [all_day_sd; sd_block_data]; % Adding each block sd in each row
        all_day_cv = [all_day_cv; cv_block_data]; % Adding each block cv in each row
        all_day_sem = [all_day_sem; sem_block_data]; % Adding each block sem in each row
        
        % Saving all trials data
        Force_Array_Name1 = sprintf('%s_Force_b%s_alltrls.txt',Folder_Name, fig1_name); % Name of the plot
        writematrix(all_force_trls,fullfile(Participant_OutputFolder_Name,Force_Array_Name1),'Delimiter','tab')
        
    end
    
    num_of_blocks = (1:length(ArrangedFileList_Force))';
    
    %% Saving all STATS variables
    avg_data_final = [num_of_blocks, all_day_avg];
    
    for dat_avg = 1:size(avg_data_final,1)
        fprintf(fid_AVG,'%.6f\t',avg_data_final(dat_avg,:));
        fprintf(fid_AVG,'\n');
    end
    fclose(fid_AVG)
   
%     sd_data_final = [num_of_blocks, all_day_sd];
%     for dat_sd = 1:size(sd_data_final,1)
%         fprintf(fid_SD,'%.6f\t',sd_data_final(dat_sd,:));
%         fprintf(fid_SD,'\n');
%     end
%     fclose(fid_SD)
%     
%     cv_data_final = [num_of_blocks, all_day_cv];
%     for dat_cv = 1:size(cv_data_final,1)
%         fprintf(fid_CV,'%.6f\t',cv_data_final(dat_cv,:));
%         fprintf(fid_CV,'\n');
%     end
%     fclose(fid_CV)
%     
%     sem_data_final = [num_of_blocks, all_day_sem];
%     for dat_sem = 1:size(sem_data_final,1)
%         fprintf(fid_SEM,'%.6f\t',sem_data_final(dat_sem,:));
%         fprintf(fid_SEM,'\n');
%     end
%     fclose(fid_SEM)
    
end

fclose('all');

toc