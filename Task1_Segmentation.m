%% BML714 LAB 4
% EMG filtering, task segmentation, and visualization
% Muscles: FDS, FCR, ECR
% Force: Dynamometer
%
% Input structure:
% Data.Trial1.EMG
% Data.Trial1.Force
% Data.Trial1.PhaseTiming
%
% Required phase timing columns:
% StartTime_s, EndTime_s, Phase
%
% Output:
% Filtered EMG data
% Segmented EMG and force data
% Separate plots for FDS, FCR, and ECR
% All results saved in GradedGrip_Results

clc;
clear;
close all;

%% 1. LOAD ORIGINAL DATA

[fileName, filePath] = uigetfile('*.mat', ...
    'Select the Graded Grip EMG Data File');

if isequal(fileName, 0)
    error('No MAT file selected. Script terminated.');
end

loadedData = load(fullfile(filePath, fileName));

if isfield(loadedData, 'Data')
    Data = loadedData.Data;
else
    error('The selected MAT file does not contain a variable named Data.');
end

disp('Original data loaded successfully.');

%% 2. SETTINGS

muscles = {'FDS', 'FCR', 'ECR'};

Fs_EMG = 1000;       % EMG sampling frequency (Hz)
Fs_Force = 10;       % Force sampling frequency (Hz)

lowCutoff = 20;      % EMG high-pass cutoff (Hz)
highCutoff = 450;    % EMG low-pass cutoff (Hz)
notchFrequency = 50; % Power-line interference (Hz)

%% 3. CREATE OUTPUT FOLDERS

outputFolder = fullfile(filePath, 'GradedGrip_Results');
plotFolder = fullfile(outputFolder, 'Plots');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

if ~exist(plotFolder, 'dir')
    mkdir(plotFolder);
end

%% 4. FILTER DESIGN

% Butterworth band-pass filter
[bBand, aBand] = butter(4, ...
    [lowCutoff highCutoff] / (Fs_EMG/2), 'bandpass');

% 50 Hz notch filter
wo = notchFrequency / (Fs_EMG/2);
bw = wo / 35;

[bNotch, aNotch] = iirnotch(wo, bw);

%% 5. IDENTIFY TRIALS

trialNames = fieldnames(Data);

trialNames = trialNames(startsWith(trialNames, 'Trial'));

if isempty(trialNames)
    error('No trial fields found. Expected Trial1, Trial2, etc.');
end

% Initialize output structures ONCE
FilteredData = struct();
SegmentedData = struct();

%% 6. PROCESS EACH TRIAL

for tr = 1:numel(trialNames)

    trialName = trialNames{tr};

    fprintf('\nProcessing %s...\n', trialName);

    %% LOAD TRIAL DATA

    EMG = Data.(trialName).EMG;
    Force = Data.(trialName).Force;
    PhaseTiming = Data.(trialName).PhaseTiming;

    % Extract time vectors
    emgTime = EMG.Time_s;
    forceTime = Force.Time_s;

    % Convert time vectors to column vectors
    emgTime = emgTime(:);
    forceTime = forceTime(:);
   

    %% 7. FILTER EMG CHANNELS

    for m = 1:numel(muscles)

        muscle = muscles{m};

        rawSignal = double(EMG.(muscle));
        rawSignal = rawSignal(:);

        % Remove DC offset
        rawSignal = rawSignal - mean(rawSignal, 'omitnan');

        % Band-pass filtering
        filteredSignal = filtfilt(bBand, aBand, rawSignal);

        % Notch filtering
        filteredSignal = filtfilt(bNotch, aNotch, filteredSignal);

        % Store filtered EMG
        FilteredData.(trialName).(muscle) = filteredSignal;

    end

    % Store EMG timestamps
    FilteredData.(trialName).EMG_Time_s = emgTime;

    %% 8. STORE FORCE DATA

    forceSignal = double(Force.Force_N);
    forceSignal = forceSignal(:);

    FilteredData.(trialName).Force_N = forceSignal;
    FilteredData.(trialName).Force_Time_s = forceTime;

    %% 9. SEGMENT EMG AND FORCE BY PHASE

    nPhases = height(PhaseTiming);

    % Preserve the original phase timing table
    SegmentedData.(trialName).PhaseTiming = PhaseTiming;

    for p = 1:nPhases

        startTime = PhaseTiming.StartTime_s(p);
        endTime = PhaseTiming.EndTime_s(p);

        phaseField = sprintf('Phase%02d', p);

        % Extract phase label
        if ismember('Phase', PhaseTiming.Properties.VariableNames)

            phaseLabel = string(PhaseTiming.Phase{p});

        elseif ismember('PhaseName', PhaseTiming.Properties.VariableNames)

            phaseLabel = string(PhaseTiming.PhaseName{p});

        else

            phaseLabel = "Phase" + string(p);

        end

        % Include the final endpoint in the last phase only
        if p < nPhases
            idxEMG = emgTime >= startTime & emgTime < endTime;
            idxForce = forceTime >= startTime & forceTime < endTime;
        else
            idxEMG = emgTime >= startTime & emgTime <= endTime;
            idxForce = forceTime >= startTime & forceTime <= endTime;
        end

        % Store phase metadata
        SegmentedData.(trialName).(phaseField).PhaseLabel = phaseLabel;

        SegmentedData.(trialName).(phaseField).StartTime_s = ...
            startTime;

        SegmentedData.(trialName).(phaseField).EndTime_s = ...
            endTime;

        % Store target percentage if available
        if ismember('TargetPercent', ...
                PhaseTiming.Properties.VariableNames)

            SegmentedData.(trialName).(phaseField).TargetPercent = ...
                PhaseTiming.TargetPercent(p);

        end

        % Store EMG timestamps
        SegmentedData.(trialName).(phaseField).EMG_Time_s = ...
            emgTime(idxEMG);

        % Store all three segmented EMG channels
        for m = 1:numel(muscles)

            muscle = muscles{m};

            filteredSignal = FilteredData.(trialName).(muscle);

            SegmentedData.(trialName).(phaseField).(muscle) = ...
                filteredSignal(idxEMG);

        end

        % Store force timestamps and values
        SegmentedData.(trialName).(phaseField).Force_Time_s = ...
            forceTime(idxForce);

        SegmentedData.(trialName).(phaseField).Force_N = ...
            forceSignal(idxForce);

    end

   %% 10. PLOT EACH EMG CHANNEL WITH FORCE

        for m = 1:numel(muscles)
        
            muscle = muscles{m};
        
            % Retrieve the signals for the current trial
            emgSignal = FilteredData.(trialName).(muscle);
            emgTimePlot = FilteredData.(trialName).EMG_Time_s;
        
            forceSignalPlot = FilteredData.(trialName).Force_N;
            forceTimePlot = FilteredData.(trialName).Force_Time_s;
        
            figure('Position', [100 100 1000 700]);
        
            %% UPPER SUBPLOT: FILTERED EMG
        
            ax1 = subplot(2,1,1);
        
            plot(ax1, emgTimePlot, emgSignal, ...
                'k', 'LineWidth', 0.8);
        
            ylabel(ax1, 'EMG Amplitude');
            title(ax1, [muscle ' - Filtered EMG']);
            grid(ax1, 'on');
            hold(ax1, 'on');
        
            %% LOWER SUBPLOT: GRIP FORCE
        
            ax2 = subplot(2,1,2);
        
            plot(ax2, forceTimePlot, forceSignalPlot, ...
                'b', 'LineWidth', 1);
        
            ylabel(ax2, 'Force (N)');
            xlabel(ax2, 'Time (s)');
            title(ax2, 'Grip Force', 'Units', 'normalized', ...
    'Position', [0.5 1.18 0]);
            grid(ax2, 'on');
            hold(ax2, 'on');
        
            %% MARK PHASE BOUNDARIES AND LABELS
        
            % Use the common time range for positioning phase labels
            xLimits = [min(emgTimePlot), max(emgTimePlot)];
            timeRange = diff(xLimits);
        
            if timeRange <= 0
                error('Invalid EMG time range in %s.', trialName);
            end
        
            for p = 1:nPhases
        
                startTime = PhaseTiming.StartTime_s(p);
                endTime = PhaseTiming.EndTime_s(p);
        
                % Extract the phase label
                if ismember('Phase', PhaseTiming.Properties.VariableNames)
        
                    if iscell(PhaseTiming.Phase)
                        phaseLabel = string(PhaseTiming.Phase{p});
                    else
                        phaseLabel = string(PhaseTiming.Phase(p));
                    end
        
                elseif ismember('PhaseName', PhaseTiming.Properties.VariableNames)
        
                    if iscell(PhaseTiming.PhaseName)
                        phaseLabel = string(PhaseTiming.PhaseName{p});
                    else
                        phaseLabel = string(PhaseTiming.PhaseName(p));
                    end
        
                else
        
                    phaseLabel = "Phase " + string(p);
        
                end
        
                % Draw phase boundaries on both subplots
                xline(ax1, startTime, '--r', 'LineWidth', 1);
                xline(ax1, endTime, '--r', 'LineWidth', 1);
        
                xline(ax2, startTime, '--r', 'LineWidth', 1);
                xline(ax2, endTime, '--r', 'LineWidth', 1);
        
                % Find the midpoint of the phase
                labelTime = (startTime + endTime) / 2;
        
                % Convert time to normalized horizontal position (0 to 1)
                xNorm = (labelTime - xLimits(1)) / timeRange;
        
                % Keep the label within the axes
                xNorm = max(0.02, min(0.98, xNorm));
        
                % Horizontal phase label above the EMG subplot
                text(ax1, xNorm, 1.04, char(phaseLabel), ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 9, ...
                    'FontWeight', 'bold', ...
                    'Clipping', 'off');
        
                % Horizontal phase label above the force subplot
                text(ax2, xNorm, 1.04, char(phaseLabel), ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 9, ...
                    'FontWeight', 'bold', ...
                    'Clipping', 'off');
        
            end
        
            %% RAISE THE UPPER SUBPLOT TITLE SLIGHTLY
        
            ax1.Title.Units = 'normalized';
            ax1.Title.Position(2) = 1.12;
        
            %% LINK X-AXES
        
            linkaxes([ax1, ax2], 'x');
        
            %% SAVE FIGURE
        
            saveas(gcf, fullfile(plotFolder, ...
                [trialName '_' muscle '_Plot.png']));
        
            savefig(gcf, fullfile(plotFolder, ...
                [trialName '_' muscle '_Plot.fig']));
        
        end
end 
%% 11. SAVE ALL FILTERED AND SEGMENTED DATA

% Save filtered EMG and original force recordings
save(fullfile(outputFolder, 'Lab4_Filtered_Data.mat'), ...
    'FilteredData', '-v7.3');

% Save all segmented EMG, force, and phase metadata
save(fullfile(outputFolder, 'Lab4_Segmented_EMG_Force_Data.mat'), ...
    'SegmentedData', '-v7.3');

%% 12. COMPLETION MESSAGE

disp(' ');
disp('Filtering, segmentation, plotting, and saving completed.');
disp(['Results saved in: ' outputFolder]);
disp(['Plots saved in: ' plotFolder]);