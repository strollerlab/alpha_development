% =============================================================================
% ALPHA BURST DEVELOPMENT PROJECT
% Script: Figures_S3_S4_S5_topomaps_fieldtrip.m
% Purpose: Creates mean scalp topographies per metric x visit (Supplementary
%          Figs. S3-S5) from the merged topological heatmap CSV, using
%          fieldtrip topoplot.
% -----------------------------------------------------------------------------
% MANUSCRIPT NOMENCLATURE Table 1 term
%   Column names are read from the merged R-pipeline CSV
%     slope                    -> Slope (aperiodic exponent)
%     offset                   -> Offset
%     peak_freq                -> Peak Frequency
%     peak_ampl                -> Peak Amplitude
%     peak_per                 -> Prop. of Peaks
%     osc_ampl                 -> Band Power
%     volt_amp / *_corrected   -> Voltage Amplitude (absolute / corrected)
%     band_amp / *_corrected   -> Band Amplitude (absolute / corrected)
%     prop_bursty_epochs       -> Prop. of Epochs w/ Burst
%     prop_bursty_cycles_burst -> Prop. of Cycles w/ Burst
%     avg_burst_duration       -> Burst Duration (consecutive burst cycles)
%     alpha_LAcH               -> Lifespan (based on Lagged Autocoherence Hilbert)
% =============================================================================

% Clear all the space and variables
clc
clear

%% Load the data
% EDIT only path2root. Output sub-folder mirrors the manuscript taxonomy (README).
% This script creates and saves the topographic heat maps for:
%   Fig S3 (aperiodic + oscillatory PSD), Fig S4 (burst properties),
%   Fig S5 (Alpha lifespan)           
%   Fig SM3 (parametrized quality metrics: R2/MAE) for Supplementary Methods. 
%   All are written here by variable name; sort the
%   R2/MAE panels into Figures/SupplementaryMethods/ when assembling SM3.
%% ==== EDIT THESE SIX LINES — nothing else in this file needs changing ====
path2data        = ''; % EDIT 
path2root        = '';  % EDIT
path2figs        = fullfile(path2root, 'SupplementaryInformation', 'GeneralInformation', 'Figures', 'IndividualTopomaps');  

electrode_layout = '~/AlphaBurstRhythm/Data/egi128_layout.sfp'; % EDIT 
path2field       = '~/toolboxes/fieldtrip-20181231/'; % EDIT
datapath         = '';   % folder holding the EEG fake dataset that hold the electrodes layout (EEG.mat)
dataname_set     = 'EEG.mat'; % name of that .mat file. For easiness, instead of providing a EEGlab dataset, here we give you a dummy MATLAB structure with the information and format necessary to run the code. 
%% ==== END OF USER CONFIGURATION ==========================================

load(dataname_set)
addpath(genpath(path2field))

EEG = pop_loadset(dataname_set, datapath); % Load the dummy dataset
EEG = pop_select(EEG, 'channel', elect);
EEG = eeglab2fieldtrip(EEG,'preprocessing');

if ~exist(fullfile(path2figs, 'individual_topomaps'), 'dir')
    mkdir(path2figs)
    mkdir(fullfile(path2figs, 'individual_topomaps'))
end 

%% Load Dataset and Prepare FieldTrip Structure

% Remove unnecessary fields from the dummy EEG
% Load Colormap
colormap_bupu = brewermap([], 'BuPu'); % Can be modified if wanted to

% --- Plot Definitions --- Create as many groups as wanted. 
groups = struct();

% Group 1: Aperiodic
groups(1).name       = 'Aperiodic_Activity';
groups(1).vars       = {'offset', 'slope'}; % Variable names in the csv
groups(1).units      = {'', ''};
groups(1).col_labels = {'Offset', 'Slope'}; % Labels that will appear on the plot
groups(1).type_activit = 'Aperiodic';

% Group 2: Oscillatory
groups(2).name       = 'Oscillatory_Activity';
groups(2).vars       = {'peak_per', 'peak_freq', 'peak_ampl', 'osc_ampl'};
groups(2).units      = {'', '', '', ''};
groups(2).col_labels = {'Prop. of Peaks', 'Peak Freq.', 'Peak Amp.', 'Band Power'};
groups(2).type_activit = 'Oscillatory';

% Group 3: Bursts 
groups(3).name       = 'Burst_Activity';
groups(3).vars       = {'volt_amp_Burst', 'volt_amp_NoBurst', 'volt_amp_corrected',...
    'band_amp_Burst', 'band_amp_NoBurst', 'band_amp_corrected',...
    'frequency_Burst', 'frequency_NoBurst',...
    'prop_bursty_cycles_burst', 'prop_bursty_epochs', 'avg_burst_duration'};
groups(3).units      = {'', '', '', '', '', '',...
    '', '', '', '', ''};
groups(3).col_labels = {'Volt. Amp.', 'Volt. Amp.', 'Corrected Volt. Amp.',...
  'Band Amp.', 'Band Amp.', 'Corrected Band Amp.',...
  'Freq.', 'Freq.',...
  'Prop. of Cycles', 'Prop. of Epochs', 'Burst Duration'};
groups(3).type_activit = 'Oscillatory_Burst';

%Group 4: LaCH
groups(4).name       = 'LAcH - Lifespan';
groups(4).vars       = {'alpha_LAcH'};
groups(4).units      = {''};
groups(4).col_labels = {'Lifespan'};
groups(4).type_activit = 'Rhythmicity';


% Group 5: Model Fit Quality   % <-- Haleigh added this
groups(5).name       = 'Model_Fit_Quality';
groups(5).vars       = {'r2value', 'mae'};
groups(5).units      = {'', ''};
groups(5).col_labels = {'R²', 'MAE'};
groups(5).type_activit = 'ModelFit';

% Conditions (Sub-Columns) - In our case it is visit age
block_keys     = unique(data.session_age); % Individual "block" values. Modify if your "Block" of interest in another one. 
block_titles   = {'1mo.', '6mo.', '12mo.', '18mo.', '30mo.', '36mo.', '42mo.', '48mo.'}; % Labels X Axis: Needs to map order of unique(data.XX)
all_vars       = data.Properties.VariableNames;

%% Main Analysis Loop
for g_idx = 1:length(groups)
    
    target_params = groups(g_idx).vars; % We locate the variable names in the analysis
    
    % --- Loop through VARIABLES inside the group (One Figure per Variable) ---
    for v_idx = 1:length(target_params)
        
        cur_param = target_params{v_idx};
        
        % Search variable name in table
        match_idx = find(contains(all_vars, cur_param));
        if isempty(match_idx), continue; end
        var_name = all_vars{match_idx(1)};
        
        % --- Figure Setup ---
        n_rows = 1;                  % Fixed to 1 row per figure: EDIT if wanted
        n_cols = length(block_keys)/n_rows; % ONLY EDIT IF THE NUMBER OF BLOCKS/AGES IS NOT DIVISIBLE BY THE NUMBER OF ROWS
        
        SUBPLOT_W = 3.5; % Individual plot size (in cm).
        SUBPLOT_H = 3.5;
        
        % Width = 8 columns + margins
        FIG_W     = SUBPLOT_W * n_cols + 4; 
        FIG_H     = SUBPLOT_H * n_rows + 2; % Much shorter height
        
        h_fig = figure('Units', 'centimeters', 'Position', [0, 0, FIG_W, FIG_H]);
        colormap(h_fig, colormap_bupu);
        h_axes = gobjects(n_rows, n_cols);
        
        % --- Calculate Common Scale for this Variable --- This is done
        % across blocks/ages to adjust the scale
        data_min = min(data.(var_name)) + 0.1* std(data.(varname)); % We give it some room so it is not saturated in case scale differences among blocks are big
        data_max = max(data.(var_name)) + 0.1* std(data.(varname));
        
        if isempty(data_min) || isempty(data_max) % THIS MUST NEVER HAPPEN
            minclim = 0; maxclim = 1;
        else
            minclim = data_min;
            maxclim = data_max;
        end

        % --- Loop through Blocks/Ages (Columns) ---
        for col_idx = 1:n_cols
            cur_age = block_keys(col_idx);
            cur_age_title = block_titles{col_idx}; % It is called cur_age_title because this study used ages. Indiferent if these are conditions, blocks, etc. 
            
            % Create Subplot (Single Row)
            h_axes(1, col_idx) = subplot(n_rows, n_cols, col_idx);
            
            % Extract Data
            data_sub_vals = data.(var_name)(data.session_age == cur_age);
            data_sub_labels = data.label(data.session_age == cur_age);
            
            dataprov = nan(length(EEG.label), 1);
            for chan_i = 1:length(EEG.label)
                row_match = find(strcmp(data_sub_labels, EEG.label{chan_i}), 1);
                if ~isempty(row_match)
                    dataprov(chan_i) = data_sub_vals(row_match);
                end
            end
            
            ft_data = EEG;
            ft_data.powspctrm = dataprov;
            
            % Topoplot Config
            cfg = [];
            cfg.marker      = 'off';
            cfg.zlim        = [minclim maxclim]; 
            cfg.ylim        = [1,1];
            cfg.layout      = electrode_layout;
            cfg.comment     = 'no';
            cfg.interactive = 'no';
            cfg.colormap    = colormap_bupu;
            cfg.colorbar    = 'no'; 
            
            ft_topoplotTFR(cfg, ft_data);
            
            % --- Labels ---
            
            % 1. Age Title (Top of every plot)
            title(cur_age_title, 'FontSize', 12, 'FontName', 'Aptos', 'FontWeight', 'bold');
            
            % 2. Variable Name (Left of first plot)
            if col_idx == 1
                text(-0.35, 0.5, groups(g_idx).col_labels{v_idx}, 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Rotation', 90, ...
                    'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Aptos');
            end
            
            % 3. Colorbar (Right of last plot)
            if col_idx == n_cols
                cbar = colorbar;
                set(cbar, 'Location', 'eastoutside');
                ylabel(cbar, groups(g_idx).units{v_idx}, 'FontSize', 8, 'FontName', 'Aptos');
                caxis([minclim maxclim]);
            end
            
            set(gca, 'FontName', 'Aptos');
        end
        
        % --- Layout Spacing: WE CROP IT SO THE WHITE SPACE BETWEEN ADJACENT PLOTS IS NOT BIG ---
        margin_l = 0.08; 
        margin_r = 0.08; 
        margin_t = 0.15; % Higher top margin for titles
        margin_b = 0.05;  
        gap_w = 0.01;     
        
        plot_width_total = 1 - margin_l - margin_r - (gap_w * (n_cols - 1));
        sp_w = plot_width_total / n_cols;
        sp_h = 1 - margin_t - margin_b; % Full height available (since only 1 row)
        
        current_y = margin_b; % Bottom up
        
        for c = 1:n_cols
            current_x = margin_l + ((c-1) * (sp_w + gap_w));
            
            if isgraphics(h_axes(1,c))
                set(h_axes(1,c), 'Position', [current_x, current_y, sp_w, sp_h]);
                
                if c == n_cols
                    cb = findobj(h_fig, 'Type', 'colorbar', 'Peer', h_axes(1,c));
                    if ~isempty(cb)
                        cb_x = current_x + sp_w + 0.01; 
                        set(cb, 'Position', [cb_x, current_y, 0.015, sp_h]);
                    end
                end
            end
        end
        
        % Save per Variable
        % Sanitizing variable name for filename
        safe_var_name = regexprep(cur_param, '[^a-zA-Z0-9]', '_');
        save_name = ['Topomap_' safe_var_name '.jpeg'];
        
        print(h_fig, fullfile(path2figs, save_name), '-djpeg', '-r600'); % Resolution of the plot. 
        close(h_fig);
        
    end % End Variable Loop
end % End Group Loop