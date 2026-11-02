%% BML714 LAB 4 - TASK 2: STATIC AND DYNAMIC HOLDING
% EMG filtering, phase segmentation, plotting, and data export
%
% EMG sampling frequency: 1000 Hz
% EMG channels: FDS, FCR, ECR
%
% Filtering:
%   20-450 Hz, 4th-order Butterworth bandpass
%   50 Hz notch filter
%
% Phases:
%   REST, RAISE, HOLD, LOWER, REST
%
% Output:
%   Static_Dynamic_Holding_Results/
%       Plots/
%       Filtered_Data/
%       Segmented_Data/
%
% NOTE:
%   This script does not use or require force data.

clear;
clc;
close all;

%% 1. SETTINGS

Fs = 1000;                  % EMG sampling frequency (Hz)
lowCutoff = 20;             % Bandpass lower cutoff (Hz)
highCutoff = 450;           % Bandpass upper cutoff (Hz)
notchFrequency = 50;        % Notch frequency (Hz)
notchQualityFactor = 35;    % Notch filter quality factor

muscleNames = {'FDS', 'FCR', 'ECR'};

outputFolder = 'Static_Dynamic_Holding_Results';
plotsFolder = fullfile(outputFolder, 'Plots');
filteredFolder = fullfile(outputFolder, 'Filtered_Data');
segmentedFolder = fullfile(outputFolder, 'Segmented_Data');

% Create output folders

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

if ~exist(plotsFolder, 'dir')
    mkdir(plotsFolder);
end

if ~exist(filteredFolder, 'dir')
    mkdir(filteredFolder);
end

if ~exist(segmentedFolder, 'dir')
    mkdir(segmentedFolder);
end

%% 2. LOAD TASK 2 DATA

[fileName, filePath] = uigetfile('*.mat', ...
    'Select the Task 2 MAT file');

if isequal(fileName, 0)
    error('No MAT file selected. Script terminated.');
end

loadedData = load(fullfile(filePath, fileName));

if isfield(loadedData, 'DataTask2')
    DataTask2 = loadedData.DataTask2;

elseif isfield(loadedData, 'Data')
    DataTask2 = loadedData.Data;

else
    error(['The selected MAT file must contain ' ...
           'DataTask2 or Data.']);
end

fprintf('\nTask 2 data loaded successfully.\n');
fprintf('File: %s\n', fileName);
fprintf('EMG sampling frequency: %d Hz\n', Fs);

%% 3. DESIGN FILTERS

% Butterworth bandpass filter: 20-450 Hz

[bBand, aBand] = butter(4, ...
    [lowCutoff highCutoff] / (Fs / 2), 'bandpass');

% 50 Hz notch filter

wo = notchFrequency / (Fs / 2);
bw = wo / notchQualityFactor;

[bNotch, aNotch] = iirnotch(wo, bw);

fprintf('\nFilters designed successfully.\n');
fprintf('Bandpass: %d-%d Hz, 4th-order Butterworth\n', ...
    lowCutoff, highCutoff);
fprintf('Notch: %d Hz\n', notchFrequency);

%% 4. IDENTIFY TRIALS

if ~isstruct(DataTask2)
    error('DataTask2 must be a structure containing trial fields.');
end

trialNames = fieldnames(DataTask2);

isTrial = startsWith(trialNames, 'Trial', ...
    'IgnoreCase', true);

trialNames = trialNames(isTrial);

if isempty(trialNames)
    error('No trial fields were found in DataTask2.');
end

fprintf('\nNumber of trials found: %d\n', numel(trialNames));

%% 5. INITIALIZE OUTPUT STRUCTURES

AllFilteredData = struct();
AllSegmentedData = struct();

%% 6. PROCESS EACH TRIAL

for trialIndex = 1:numel(trialNames)

    trialName = trialNames{trialIndex};

    fprintf('\n========================================\n');
    fprintf('Processing %s...\n', trialName);
    fprintf('========================================\n');

    trialData = DataTask2.(trialName);

    %% 6.1 CHECK REQUIRED FIELDS

    % Only EMG and PhaseTiming are required.
    % Force data are not accessed anywhere in this script.

    if ~isfield(trialData, 'EMG')
        warning('%s does not contain EMG data. Skipping.', trialName);
        continue;
    end

    if ~isfield(trialData, 'PhaseTiming')
        warning('%s does not contain PhaseTiming. Skipping.', trialName);
        continue;
    end

    EMG = trialData.EMG;
    PhaseTiming = trialData.PhaseTiming;

    if ~istable(EMG)
        warning('%s EMG data are not stored as a table. Skipping.', ...
            trialName);
        continue;
    end

    if ~istable(PhaseTiming)
        warning('%s PhaseTiming is not a table. Skipping.', trialName);
        continue;
    end

    %% 6.2 VALIDATE EMG COLUMNS

    requiredEMGColumns = [{'Time_s'}, muscleNames];

    if ~all(ismember(requiredEMGColumns, ...
            EMG.Properties.VariableNames))

        error('%s is missing one or more required EMG columns.', ...
            trialName);
    end

    originalEMGTime = double(EMG.Time_s(:));

    nSamples = height(EMG);

    if nSamples < 10
        warning('%s has too few EMG samples. Skipping.', trialName);
        continue;
    end

    % Reconstruct the time vector using the known 1000 Hz
    % acquisition rate.
    %
    % Assumption: EMG samples were acquired uniformly at 1000 Hz.
    % Original timestamps are retained separately for inspection.

    emgTime = (0:nSamples-1)' / Fs;

    fprintf('Number of EMG samples: %d\n', nSamples);
    fprintf('Expected sampling frequency: %d Hz\n', Fs);
    fprintf('Reconstructed EMG duration: %.3f s\n', emgTime(end));

    % Display original timestamp diagnostics.

    originalDt = diff(originalEMGTime);
    finiteDt = originalDt(isfinite(originalDt));

    if ~isempty(finiteDt)

        fprintf('Original timestamp median difference: %.9f s\n', ...
            median(finiteDt));

        fprintf('Original timestamp differences <= 0: %d\n', ...
            sum(finiteDt <= 0));

    end

    %% 6.3 FILTER EMG CHANNELS

    filteredEMG = struct();

    for muscleIndex = 1:numel(muscleNames)

        muscle = muscleNames{muscleIndex};

        rawSignal = double(EMG.(muscle)(:));

        % Identify finite samples.

        validSamples = isfinite(rawSignal);

        if sum(validSamples) < 10

            warning('Insufficient valid samples for %s in %s.', ...
                muscle, trialName);

            continue;
        end

        % Replace non-finite samples by linear interpolation,
        % if necessary.

        if any(~validSamples)

            sampleIndex = (1:nSamples)';

            rawSignal(~validSamples) = interp1( ...
                sampleIndex(validSamples), ...
                rawSignal(validSamples), ...
                sampleIndex(~validSamples), ...
                'linear', 'extrap');

            fprintf('%s: replaced %d non-finite samples.\n', ...
                muscle, sum(~validSamples));

        end

        % Verify signal validity.

        if any(~isfinite(rawSignal))
            error('Non-finite samples remain in %s.', muscle);
        end

        % Remove the DC offset.

        rawSignal = rawSignal - mean(rawSignal);

        % Apply the Butterworth bandpass filter.

        try

            bandpassSignal = filtfilt( ...
                bBand, aBand, rawSignal);

        catch ME

            error(['Bandpass filtering failed for %s in %s: %s\n' ...
                   'Stop processing and check filter stability.'], ...
                muscle, trialName, ME.message);

        end

        % Check bandpass output.

        if any(~isfinite(bandpassSignal))

            error(['Bandpass filtering produced NaN/Inf values ' ...
                   'for %s in %s.'], muscle, trialName);

        end

        % Apply the 50 Hz notch filter.

        try

            finalSignal = filtfilt( ...
                bNotch, aNotch, bandpassSignal);

        catch ME

            error('Notch filtering failed for %s in %s: %s', ...
                muscle, trialName, ME.message);

        end

        % Check final filtered signal.

        if any(~isfinite(finalSignal))

            error(['Notch filtering produced NaN/Inf values ' ...
                   'for %s in %s.'], muscle, trialName);

        end

        % Store filtered signal.

        filteredEMG.(muscle) = finalSignal;

        % Display signal diagnostics.

        fprintf(['%s | %s | Samples: %d | Min: %.5f | ' ...
                 'Max: %.5f | SD: %.5f\n'], ...
            trialName, muscle, numel(finalSignal), ...
            min(finalSignal), max(finalSignal), std(finalSignal));

    end

    % Ensure all three muscles were processed.

    if ~all(isfield(filteredEMG, muscleNames))

        warning(['One or more muscle channels could not be filtered ' ...
                 'in %s. Skipping segmentation and plotting.'], ...
                 trialName);

        continue;

    end

    %% 6.4 VALIDATE PHASE TIMING

    requiredPhaseColumns = ...
        {'StartTime_s', 'EndTime_s', 'Phase'};

    if ~all(ismember(requiredPhaseColumns, ...
            PhaseTiming.Properties.VariableNames))

        error('%s is missing required phase timing columns.', ...
            trialName);

    end

    phaseStart = double(PhaseTiming.StartTime_s(:));
    phaseEnd = double(PhaseTiming.EndTime_s(:));
    phaseNames = string(PhaseTiming.Phase(:));

    if any(~isfinite(phaseStart)) || any(~isfinite(phaseEnd))

        error('%s contains non-finite phase boundaries.', ...
            trialName);

    end

    if any(phaseEnd <= phaseStart)

        error('%s contains invalid phase intervals.', ...
            trialName);

    end

    % Sort phases chronologically.

    [phaseStart, phaseOrder] = sort(phaseStart);

    phaseEnd = phaseEnd(phaseOrder);
    phaseNames = phaseNames(phaseOrder);

    % Warn if phases fall outside the reconstructed time range.

    if any(phaseStart < emgTime(1)) || ...
            any(phaseEnd > emgTime(end))

        warning(['Some phase boundaries in %s extend beyond ' ...
                 'the reconstructed EMG time range. Verify the ' ...
                 'time origin and phase timing.'], trialName);

    end

    %% 6.5 SEGMENT EMG BY PHASE

    segmentedTrial = struct();

    for phaseIndex = 1:numel(phaseNames)

        currentPhase = char(phaseNames(phaseIndex));

        startTime = phaseStart(phaseIndex);
        endTime = phaseEnd(phaseIndex);

        % Identify samples belonging to this phase.

        phaseMask = emgTime >= startTime & ...
                    emgTime <= endTime;

        % Make a valid field name and preserve repeated phases.
        % Example: REST_01 and REST_05.

        phaseField = matlab.lang.makeValidName(currentPhase);

        phaseField = sprintf('%s_%02d', ...
            phaseField, phaseIndex);

        % Store phase metadata.

        segmentedTrial.(phaseField).Phase = currentPhase;

        segmentedTrial.(phaseField).StartTime_s = startTime;

        segmentedTrial.(phaseField).EndTime_s = endTime;

        segmentedTrial.(phaseField).EMG_Time_s = ...
            emgTime(phaseMask);

        % Store the segmented EMG for each muscle.

        for muscleIndex = 1:numel(muscleNames)

            muscle = muscleNames{muscleIndex};

            signal = filteredEMG.(muscle);

            segmentedTrial.(phaseField).(muscle) = ...
                signal(phaseMask);

        end

        fprintf('%s | %s | %.2f-%.2f s | %d samples\n', ...
            trialName, currentPhase, startTime, endTime, ...
            sum(phaseMask));

    end

    %% 6.6 STORE FULL FILTERED DATA

    AllFilteredData.(trialName).EMG_Time_s = emgTime;

    AllFilteredData.(trialName).Original_EMG_Time_s = ...
        originalEMGTime;

    AllFilteredData.(trialName).Fs = Fs;

    AllFilteredData.(trialName).PhaseTiming = PhaseTiming;

    for muscleIndex = 1:numel(muscleNames)

        muscle = muscleNames{muscleIndex};

        AllFilteredData.(trialName).(muscle) = ...
            filteredEMG.(muscle);

    end

    %% 6.7 STORE SEGMENTED DATA

    AllSegmentedData.(trialName) = segmentedTrial;

    %% 6.8 SAVE FILTERED DATA FOR THIS TRIAL

    save(fullfile(filteredFolder, ...
        [trialName '_Filtered_Data.mat']), ...
        'filteredEMG', 'emgTime', 'originalEMGTime', ...
        'PhaseTiming', 'Fs');

    %% 6.9 SAVE SEGMENTED DATA FOR THIS TRIAL

    save(fullfile(segmentedFolder, ...
        [trialName '_Segmented_Data.mat']), ...
        'segmentedTrial', 'Fs');

    %% 6.11 PLOT FILTERED EMG WITH PHASE LABELS

    muscles = {'FDS', 'FCR', 'ECR'};
    yLabels = {'FDS EMG', 'FCR EMG', 'ECR EMG'};
    
    fig = figure( ...
        'Name', ['Task 2 - ' trialName], ...
        'Color', 'w', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.03 0.90 0.94]);
    
    tl = tiledlayout(fig, 3, 1, ...
        'TileSpacing', 'loose', ...
        'Padding', 'loose');
    
    for k = 1:3
    
        ax = nexttile(tl);
    
        % Plot filtered EMG
        plot(ax, emgTime, filteredEMG.(muscles{k}), ...
            'Color', [0.15 0.15 0.15], ...
            'LineWidth', 0.45);
    
        hold(ax, 'on');
    
        % Draw phase start boundaries
        for p = 1:numel(phaseStart)
    
            boundary = phaseStart(p);
    
            if boundary >= emgTime(1) && boundary <= emgTime(end)
                xline(ax, boundary, '--', ...
                    'Color', [0.85 0.20 0.20], ...
                    'LineWidth', 1.0, ...
                    'HandleVisibility', 'off');
            end
        end
    
        % Draw final phase end boundary
        finalBoundary = max(phaseEnd);
    
        if finalBoundary >= emgTime(1) && finalBoundary <= emgTime(end)
            xline(ax, finalBoundary, '--', ...
                'Color', [0.85 0.20 0.20], ...
                'LineWidth', 1.0, ...
                'HandleVisibility', 'off');
        end
    
        % Add phase labels centred above each phase
        yLimits = ylim(ax);
        labelY = yLimits(2) + 0.07 * (yLimits(2) - yLimits(1));
    
        for p = 1:numel(phaseNames)
    
            phaseMidpoint = (phaseStart(p) + phaseEnd(p)) / 2;
    
            if phaseMidpoint >= emgTime(1) && phaseMidpoint <= emgTime(end)
    
                text(ax, phaseMidpoint, labelY, char(phaseNames(p)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 9, ...
                    'FontWeight', 'bold', ...
                    'Color', [0.25 0.25 0.25], ...
                    'Clipping', 'off');
    
            end
        end
    
        % Axis formatting
        ylabel(ax, yLabels{k}, ...
            'FontSize', 11, ...
            'FontWeight', 'bold');
    
        xlim(ax, [emgTime(1), emgTime(end)]);
    
        grid(ax, 'on');
        box(ax, 'on');
    
        ax.FontSize = 10;
        ax.LineWidth = 0.8;
        ax.GridAlpha = 0.15;
        ax.Layer = 'bottom';
        ax.TickDir = 'out';
    
        % Remove subplot titles
        title(ax, '');
    
        % Hide time tick labels on upper plots
        if k < 3
            ax.XTickLabel = [];
        end
    
        hold(ax, 'off');
    
    end
    
    % Shared x-axis label
    xlabel(tl, 'Time (s)', ...
        'FontSize', 11, ...
        'FontWeight', 'bold');
    
    % No overall title
    title(tl, '');
    
    %% 6.12 SAVE PLOTS
    
    if ~exist(plotsFolder, 'dir')
        mkdir(plotsFolder);
    end
    
    exportgraphics(fig, ...
        fullfile(plotsFolder, [trialName '_EMG_Plot.png']), ...
        'Resolution', 300);
    
    savefig(fig, ...
        fullfile(plotsFolder, [trialName '_EMG_Plot.fig']));
    
    close(fig);
    
    fprintf('Saved EMG plot with phase labels for %s.\n', trialName);
end

%% 7. SAVE COMBINED OUTPUT FILES

save(fullfile(outputFolder, ...
    'All_Filtered_Static_Dynamic_Holding_Data.mat'), ...
    'AllFilteredData', 'Fs', 'muscleNames', ...
    'lowCutoff', 'highCutoff', 'notchFrequency');

save(fullfile(outputFolder, ...
    'All_Segmented_Static_Dynamic_Holding_Data.mat'), ...
    'AllSegmentedData', 'Fs', 'muscleNames');

%% 8. COMPLETION MESSAGE

fprintf('\n========================================\n');
fprintf('TASK 2 EMG PROCESSING COMPLETE\n');
fprintf('========================================\n');

fprintf('\nResults folder:\n%s\n', outputFolder);

fprintf('\nOutput subfolders:\n');
fprintf('1. Plots\n');
fprintf('2. Filtered_Data\n');
fprintf('3. Segmented_Data\n');

fprintf('\nCombined output files:\n');
fprintf('1. All_Filtered_Static_Dynamic_Holding_Data.mat\n');
fprintf('2. All_Segmented_Static_Dynamic_Holding_Data.mat\n');

fprintf('\nRecommended file for phase-specific analysis:\n');
fprintf('All_Segmented_Static_Dynamic_Holding_Data.mat\n');

fprintf('\nRecommended file for continuous EMG analysis:\n');
fprintf('All_Filtered_Static_Dynamic_Holding_Data.mat\n');