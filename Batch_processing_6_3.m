%% Compute SSM data with and without NDVI stratification; SSM negative anomalies as an example
clear
clc
% Define file paths
lst_filepath = "F:\GLDAS_SSM\";  % SSM filepath; LST NDVI
ndvi_filepath = "F:\NDVI_mat\"; % NDVI filepath
output_folder = 'F:\Output_tables_figs';
eq_csvdata = "F:\EQs_anomaly_plot\EQs_query_screened_with_land_info.csv";
eq_query = readtable(eq_csvdata);

% Define radius ranges
radius_ranges = {
    '0-200km Circles', 0, 200;
    '200-300km Rings', 200, 300;
    '300-400km Rings', 300, 400;
    '400-500km Rings', 400, 500;
    '0-1000km Circles', 0, 1000;
    '0-1200km Circles', 0, 1200 
};

% Initialize grid and distance calculations
resolution = 0.05; 
rows = 3600;
cols = 7200;
R = 6371; 
lat_grid = linspace(90 - resolution/2, -90 + resolution/2, rows);
lon_grid = linspace(-180 + resolution/2, 180 - resolution/2, cols);
[LON, LAT] = meshgrid(lon_grid, lat_grid);

% Process loop 
for ij= 1:height(eq_query)
     
    % Parameter setup and time range definition
    EQ_time=eq_query.time(ij);
    EQ_date = strrep(EQ_time{1}(1:10), '-', '');
    EQ_year = str2num(EQ_date(1:4));
    EQ_month = str2num(EQ_date(5:6));
    EQ_lat = eq_query.latitude(ij); 
    EQ_lon = eq_query.longitude(ij);
    EQ_depth = eq_query.depth(ij);
    EQ_mag = eq_query.mag(ij);
    [dataset_start, dataset_end] = set_dataset_dates(EQ_date); % STL start and end input
    
    % Read R1_Land_Prop value
    land_prop = eq_query.R1_Land_Prop(ij);
    
    % Check condition: if land proportion < 40, skip; empirical threshold to retain continental coastal earthquakes (e.g., Chile)
    if land_prop < 40 || EQ_lon > 157
        continue; % Skip to next iteration
    end
    fprintf('==== Processing earthquake data %d/%d, R1 Land Proportion: %d%% ====\n', ij, height(eq_query), land_prop);
    
    baseline_length = 48; % 48-month baseline period, 24-month impact period
    baseline_start = datetime(EQ_year, EQ_month, 1) - calmonths(baseline_length + 12);
    baseline_end = dateshift(datetime(EQ_year, EQ_month, 1) - calmonths(13), 'end', 'month');
    impact_start = datetime(EQ_year, EQ_month, 1) - calmonths(12);
    impact_end = dateshift(datetime(EQ_year, EQ_month, 1) + calmonths(11), 'end', 'month');
            
    dlat = deg2rad(LAT - EQ_lat);
    dlon = deg2rad(LON - EQ_lon);
    a = sin(dlat/2).^2 + cos(deg2rad(EQ_lat)) * cos(deg2rad(LAT)) .* sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distances = R * c; 
    
    % Create 1200km mask to define calculation extent
    mask_1200km = distances <= 1200;
    
    % Find minimum bounding box containing the 1200km range
    [row_indices, col_indices] = find(mask_1200km);
    if isempty(row_indices)
        fprintf('Warning: No valid data within 1200km of earthquake location %f, %f. Skipping process.\n', EQ_lat, EQ_lon);
        continue;
    end
    
    min_row = min(row_indices);
    max_row = max(row_indices);
    min_col = min(col_indices);
    max_col = max(col_indices);
    
    % Extract sub-region dimensions
    sub_rows = max_row - min_row + 1;
    sub_cols = max_col - min_col + 1;
    
    % Create distance masks (using sub-region dimensions)
    masks = cell(size(radius_ranges, 1), 1);
    sub_distances = distances(min_row:max_row, min_col:max_col);
    
    for i = 1:size(radius_ranges, 1)
        min_radius = radius_ranges{i, 2};
        max_radius = radius_ranges{i, 3};
        
        if min_radius == 0
            masks{i} = sub_distances <= max_radius;
        else
            masks{i} = sub_distances > min_radius & sub_distances <= max_radius;
        end
    end
    
    % Retrieve all NDVI .mat files to be processed
    all_files = dir(fullfile(ndvi_filepath, 'NDVI_*.mat'));
    all_dates = [];
    
    for i = 1:length(all_files)
        filename = all_files(i).name;
        % Extract date information from filename: NDVI_YYYY_DDD.mat
        date_parts = strsplit(filename(6:end-4), '_'); % Remove "NDVI_" and ".mat"
        if length(date_parts) == 2
            yeaR = str2double(date_parts{1});
            DOY = str2double(date_parts{2});
            date_val = datetime(yeaR, 1, 1) + days(DOY - 1);
            all_dates = [all_dates; date_val];
        end
    end
    
    % Compute NDVI temporal mask
    % fprintf('Starting calculation of NDVI temporal mask...\n');
    ndvi_mask_start = datetime(EQ_year, EQ_month, 1) - calmonths(12);
    ndvi_mask_end = datetime(EQ_year, EQ_month, 1) + calmonths(12) - days(1);
    
    ndvi_mask_indices = (all_dates >= ndvi_mask_start) & (all_dates <= ndvi_mask_end);
    ndvi_mask_files = all_files(ndvi_mask_indices);
    ndvi_mask_dates = all_dates(ndvi_mask_indices);
    
    % Initialize NDVI accumulation arrays (using sub-region dimensions)
    ndvi_accumulated = zeros(sub_rows, sub_cols);
    ndvi_count = zeros(sub_rows, sub_cols);
    fprintf('Processing NDVI files %s ~ %s\n',ndvi_mask_start, ndvi_mask_end);
    for i = 1:length(ndvi_mask_files)
        filename = fullfile(ndvi_filepath, ndvi_mask_files(i).name);
        % fprintf('Processing NDVI file: %s\n', ndvi_mask_files(i).name);
        
        try
            % Load NDVI data from .mat file
            data = load(filename);
            NDVI_data = data.processed_data; 
            
            % Extract sub-region data
            NDVI_sub = NDVI_data(min_row:max_row, min_col:max_col);
                   
            spatial_mask = sub_distances <= 1000;
            valid_mask = spatial_mask & NDVI_sub ~= 0 & ~isnan(NDVI_sub) & NDVI_sub > -1;
            
            ndvi_accumulated(valid_mask) = ndvi_accumulated(valid_mask) + NDVI_sub(valid_mask);
            ndvi_count(valid_mask) = ndvi_count(valid_mask) + 1;
            
        catch ME
            fprintf('Failed to process NDVI file: %s, Error: %s\n', filename, ME.message);
            continue;
        end
    end
    
    % Set NDVI ranges - Adding case for no NDVI distinction
    NDVI_ranges = [0, 1.0;  % New: No NDVI distinction, covering full range
                   0, 0.2; 
                   0.2, 0.4; 
                   0.4, 0.6; 
                   0.6, 1.0];
    
    % Initialize result storage - Now 5 rows (1 without NDVI + 4 with NDVI)
    final_results = cell(8, 1); % 8 types of standardized anomaly indices
    for i = 1:8
        final_results{i} = zeros(5, 72); % 5 rows (1 no NDVI + 4 with NDVI) x 72 columns (4 radii x 3 intervals x 6 metrics)
    end
    
    % Pre-load all SSM data into a 3D matrix (using sub-region dimensions)
    fprintf('Reading SSM data %s ~ %s\n', dataset_start, dataset_end);
    
    % Retrieve all SSM .mat files to be processed
    lst_all_files = dir(fullfile(lst_filepath, 'SSM_*.mat'));
    lst_all_dates = [];
    
    for i = 1:length(lst_all_files)
        filename = lst_all_files(i).name;
        % Extract date information from filename: SSM_YYYY_DDD.mat
        date_parts = strsplit(filename(5:end-4), '_'); 
        if length(date_parts) == 2
            yeaR = str2double(date_parts{1});
            DOY = str2double(date_parts{2});
            date_val = datetime(yeaR, 1, 1) + days(DOY - 1);
            lst_all_dates = [lst_all_dates; date_val]; % SSM
        end
    end
    valid_indices = (lst_all_dates >= dataset_start) & (lst_all_dates <= dataset_end);
    process_dates = lst_all_dates(valid_indices);
    process_files = lst_all_files(valid_indices);
    
    % Initialize 3D matrix to store SSM data (using sub-region dimensions)
    num_days = length(process_dates);
    SSM_data_3D = zeros(sub_rows, sub_cols, num_days);
    
    % Read all SSM data (extracting only sub-region)
    for file_idx = 1:num_days
        filename = fullfile(lst_filepath, process_files(file_idx).name);
        
        try
            data = load(filename);
            % Extract data for sub-region only
            SSM_data_3D(:, :, file_idx) = data.SSM(min_row:max_row, min_col:max_col);
        catch ME
            fprintf('Failed to read SSM file: %s, Error: %s\n', filename, ME.message);
            SSM_data_3D(:, :, file_idx) = NaN;
        end
    end
    
    % Loop through 5 NDVI ranges (1st is the no-NDVI distinction case)
    for ndvi_range_idx = 1:5
        if ndvi_range_idx == 1
            fprintf('== Processing no NDVI distinction range (1/5) ==\n');
        else
            fprintf('== Processing NDVI range %.1f-%.1f (%d/5) ==\n', ...
                NDVI_ranges(ndvi_range_idx, 1), NDVI_ranges(ndvi_range_idx, 2), ndvi_range_idx);
        end
        
        NDVI_threshold_low = NDVI_ranges(ndvi_range_idx, 1);
        NDVI_threshold_high = NDVI_ranges(ndvi_range_idx, 2);
    
        % Calculate NDVI mean
        ndvi_mean = zeros(sub_rows, sub_cols);
        valid_pixels = ndvi_count > 0;
        ndvi_mean(valid_pixels) = ndvi_accumulated(valid_pixels) ./ ndvi_count(valid_pixels);
    
        % Create NDVI threshold mask
        if ndvi_range_idx == 1
            % No NDVI distinction: use all valid pixels
            ndvi_threshold_mask = valid_pixels;
        else
            % Distinguish NDVI: use threshold ranges
            ndvi_threshold_mask = (ndvi_mean >= NDVI_threshold_low) & (ndvi_mean <= NDVI_threshold_high) & valid_pixels;
        end
    
        % Initialize result structure
        results = struct();
        results.daily_means = [];
        results.daily_stds = [];
        results.daily_counts = [];
        results.dates = [];
        
        % Process each SSM file (retrieve data directly from 3D matrix)
        fprintf('Processing SSM data %s ~ %s\n', dataset_start, dataset_end);
        for file_idx = 1:num_days
            current_date = process_dates(file_idx);
            
            % Retrieve SSM data directly from 3D matrix
            SSM_data = SSM_data_3D(:, :, file_idx);
            
            daily_means = zeros(1, 5);
            daily_stds = zeros(1, 5);
            daily_counts = zeros(1, 5);
            
            for radius_idx = 1:5 % Process only the first 5 radius ranges (0-200, 200-300, 300-400, 400-500, 0-1000)
                mask = masks{radius_idx} & ndvi_threshold_mask;
                valid_mask = mask & SSM_data > 0 & ~isnan(SSM_data);
                valid_data = SSM_data(valid_mask);
                
                if ~isempty(valid_data)
                    daily_means(radius_idx) = mean(valid_data);
                    daily_stds(radius_idx) = std(valid_data);
                    daily_counts(radius_idx) = numel(valid_data);
                else
                    daily_means(radius_idx) = NaN;
                    daily_stds(radius_idx) = NaN;
                    daily_counts(radius_idx) = 0;
                end
            end
            
            results.daily_means = [results.daily_means; daily_means];
            results.daily_stds = [results.daily_stds; daily_stds];
            results.daily_counts = [results.daily_counts; daily_counts];
            results.dates = [results.dates; current_date];
        end
    
        % Calculate baseline monthly statistics (Raw SSM)
        baseline_mask = (results.dates >= baseline_start) & (results.dates <= baseline_end);
        baseline_dates = results.dates(baseline_mask);
        baseline_means = results.daily_means(baseline_mask, :);
        baseline_months = month(baseline_dates);
    
        baseline_monthly_means = zeros(12, 5);
        baseline_monthly_stds = zeros(12, 5);
        
        for m = 1:12
            month_mask = baseline_months == m;
            for radius_idx = 1:5
                month_data = baseline_means(month_mask, radius_idx);
                valid_data = month_data(~isnan(month_data));
                
                if ~isempty(valid_data)
                    baseline_monthly_means(m, radius_idx) = mean(valid_data);
                    baseline_monthly_stds(m, radius_idx) = std(valid_data);
                else
                    baseline_monthly_means(m, radius_idx) = NaN;
                    baseline_monthly_stds(m, radius_idx) = NaN;
                end
            end
        end
    
        % Calculate impact period monthly statistics (Raw SSM)
        impact_mask = (results.dates >= impact_start) & (results.dates <= impact_end);
        impact_dates = results.dates(impact_mask);
        impact_means = results.daily_means(impact_mask, :);
        impact_years = year(impact_dates);
        impact_months = month(impact_dates);
        
        unique_month_pairs = unique([impact_years, impact_months], 'rows');
        
        % Initialize 4 types of standardized anomaly indices (modified to negative anomalies)
        anomalies1 = []; % Type 1: Standardized anomaly of raw monthly mean
        anomalies2 = []; % Type 2: Standardized anomaly of monthly mean of daily means below baseline mean
        anomalies3 = []; % Type 3: Standardized anomaly of STL residual monthly mean
        anomalies4 = []; % Type 4: Standardized anomaly of STL negative residual monthly mean
        
        impact_month_offsets = [];
        
        for i = 1:size(unique_month_pairs, 1)
            year_val = unique_month_pairs(i, 1);
            month_val = unique_month_pairs(i, 2);
            
            month_mask = (impact_years == year_val) & (impact_months == month_val);
            
            monthly_mean1 = zeros(1, 5); % Raw monthly mean
            monthly_mean2 = zeros(1, 5); % Monthly mean of daily values below baseline mean
            monthly_anomaly1 = zeros(1, 5); % Type 1 standardized anomaly
            monthly_anomaly2 = zeros(1, 5); % Type 2 standardized anomaly
            
            for radius_idx = 1:5
                month_data = impact_means(month_mask, radius_idx);
                valid_data = month_data(~isnan(month_data));
                
                if ~isempty(valid_data)
                    % Raw monthly mean
                    monthly_mean1(radius_idx) = mean(valid_data);
                    
                    % Monthly mean of daily values below baseline mean
                    baseline_mean_val = baseline_monthly_means(month_val, radius_idx);
                    below_mean_data = valid_data(valid_data < baseline_mean_val);
                    
                    if ~isempty(below_mean_data)
                        monthly_mean2(radius_idx) = mean(below_mean_data);
                    else
                        monthly_mean2(radius_idx) = NaN;
                    end
                    
                    % Calculate standardized anomaly
                    baseline_mean = baseline_monthly_means(month_val, radius_idx);
                    baseline_std = baseline_monthly_stds(month_val, radius_idx);
                    
                    if baseline_std > 0 && ~isnan(baseline_mean)
                        if ~isnan(monthly_mean1(radius_idx))
                            monthly_anomaly1(radius_idx) = (monthly_mean1(radius_idx) - baseline_mean) / baseline_std;
                        end
                        
                        if ~isnan(monthly_mean2(radius_idx))
                            monthly_anomaly2(radius_idx) = (monthly_mean2(radius_idx) - baseline_mean) / baseline_std;
                        end
                    end
                end
            end
            
            % Calculate month offset
            total_months = (year_val - EQ_year) * 12 + (month_val - EQ_month); % month offset
            month_offset = total_months;
            
            impact_month_offsets = [impact_month_offsets; month_offset];
            anomalies1 = [anomalies1; monthly_anomaly1];
            anomalies2 = [anomalies2; monthly_anomaly2];
        end
    
        % STL Decomposition: Decompose daily time series for each radius range
        fprintf('Executing STL decomposition...\n');
        results.daily_residuals = zeros(size(results.daily_means));
        
        for radius_idx = 1:5
            daily_data = results.daily_means(:, radius_idx);
            dates = results.dates;
            
            tt = timetable(dates, daily_data);
            tt.Properties.VariableNames = {'LST'};
            tt_filled = fillmissing(tt, 'linear');
            
            try
                stl_result = stl(tt_filled.LST, 'Period', 365, 'Seasonal', 13, 'Trend', 549, 'Robust', true);
                results.daily_residuals(:, radius_idx) = stl_result.residual;
            catch
                results.daily_residuals(:, radius_idx) = NaN;
            end
        end
    
        % Calculate baseline monthly statistics (STL residuals)
        baseline_residuals = results.daily_residuals(baseline_mask, :);
        
        baseline_residual_means = zeros(12, 5);
        baseline_residual_stds = zeros(12, 5);
        
        for m = 1:12
            month_mask = baseline_months == m;
            for radius_idx = 1:5
                month_data = baseline_residuals(month_mask, radius_idx);
                valid_data = month_data(~isnan(month_data));
                
                if ~isempty(valid_data)
                    baseline_residual_means(m, radius_idx) = mean(valid_data);
                    baseline_residual_stds(m, radius_idx) = std(valid_data);
                else
                    baseline_residual_means(m, radius_idx) = NaN;
                    baseline_residual_stds(m, radius_idx) = NaN;
                end
            end
        end
    
        % Calculate impact period monthly statistics (STL residuals)
        impact_residuals = results.daily_residuals(impact_mask, :);
        
        % Compute Type 3 and Type 4 standardized anomaly indices
        for i = 1:size(unique_month_pairs, 1)
            year_val = unique_month_pairs(i, 1);
            month_val = unique_month_pairs(i, 2);
            
            month_mask = (impact_years == year_val) & (impact_months == month_val);
            
            monthly_mean3 = zeros(1, 5); % STL residual monthly mean
            monthly_mean4 = zeros(1, 5); % STL negative residual monthly mean
            monthly_anomaly3 = zeros(1, 5); % Type 3 standardized anomaly
            monthly_anomaly4 = zeros(1, 5); % Type 4 standardized anomaly
            
            for radius_idx = 1:5
                month_data = impact_residuals(month_mask, radius_idx);
                valid_data = month_data(~isnan(month_data));
                
                if ~isempty(valid_data)
                    % STL residual monthly mean
                    monthly_mean3(radius_idx) = mean(valid_data);
                    
                    % STL negative residual monthly mean
                    baseline_residual_mean_val = baseline_residual_means(month_val, radius_idx);
                    below_residual_data = valid_data(valid_data < baseline_residual_mean_val);
                    
                    if ~isempty(below_residual_data)
                        monthly_mean4(radius_idx) = mean(below_residual_data);
                    else
                        monthly_mean4(radius_idx) = NaN;
                    end
                    
                    % Calculate standardized anomaly
                    baseline_mean = baseline_residual_means(month_val, radius_idx);
                    baseline_std = baseline_residual_stds(month_val, radius_idx);
                    
                    if baseline_std > 0 && ~isnan(baseline_mean)
                        if ~isnan(monthly_mean3(radius_idx))
                            monthly_anomaly3(radius_idx) = (monthly_mean3(radius_idx) - baseline_mean) / baseline_std;
                        end
                        
                        if ~isnan(monthly_mean4(radius_idx))
                            monthly_anomaly4(radius_idx) = (monthly_mean4(radius_idx) - baseline_mean) / baseline_std;
                        end
                    end
                end
            end
            
            anomalies3 = [anomalies3; monthly_anomaly3];
            anomalies4 = [anomalies4; monthly_anomaly4];
        end
    
        % Sort by month offset
        [month_offset_sorted, sort_idx] = sort(impact_month_offsets);
        anomalies1_sorted = anomalies1(sort_idx, :);
        anomalies2_sorted = anomalies2(sort_idx, :);
        anomalies3_sorted = anomalies3(sort_idx, :);
        anomalies4_sorted = anomalies4(sort_idx, :);
    
        % Function to calculate 6 statistical metrics (modified for negative anomaly statistics)
        calculateStats = @(anomalies, mask) struct(...
            'MeanNegativeAnomaly', mean(anomalies(anomalies < 0 & mask)), ...
            'CumulativeNegativeIntensity', sum(anomalies(anomalies < 0 & mask)), ...
            'ProportionSignificant', sum(anomalies < 0 & mask) / sum(mask), ...
            'MaxNegativeAnomaly', min([anomalies(anomalies < 0 & mask); 0]), ...
            'MeanAnomaly', mean(anomalies(mask)), ...
            'MaxNegativeSequence', maxConsecutiveNegative(anomalies, mask) ...
        );
        
        % Calculate statistical metrics for the 4 anomaly indices
        all_anomalies = {anomalies1_sorted, anomalies2_sorted, anomalies3_sorted, anomalies4_sorted};
        
        for anomaly_idx = 1:4
            current_anomalies = all_anomalies{anomaly_idx};
            stats_matrix = zeros(4, 18); % 4 radii x 18 metrics (3 intervals x 6 metrics)
            
            for radius_idx = 1:4
                % Define masks for three time periods
                pre_mask = (month_offset_sorted >= -11 & month_offset_sorted <= 0); % 12 months pre-earthquake
                post_mask = (month_offset_sorted >= 0 & month_offset_sorted <= 11);  % 12 months post-earthquake
                full_mask = (month_offset_sorted >= -11 & month_offset_sorted <= 11); % Full 24-month period
                
                % Calculate statistics for the three time periods
                pre_stats = calculateStats(current_anomalies(:, radius_idx), pre_mask);
                post_stats = calculateStats(current_anomalies(:, radius_idx), post_mask);
                full_stats = calculateStats(current_anomalies(:, radius_idx), full_mask);
                
                % Store in matrix
                stats_matrix(radius_idx, 1:6) = struct2array(pre_stats);
                stats_matrix(radius_idx, 7:12) = struct2array(post_stats);
                stats_matrix(radius_idx, 13:18) = struct2array(full_stats);
            end
            
            % Store in final results
            final_results{anomaly_idx}(ndvi_range_idx, :) = reshape(stats_matrix', 1, []);
        end
    
        % Two-Way Fixed Effects Model Analysis (Anomaly Indices 5-8)
        for anomaly_idx = 1:4
            current_anomalies = all_anomalies{anomaly_idx};
            stats_matrix = zeros(4, 18); % 4 radii x 18 metrics (3 intervals x 6 metrics)
            
            for radius_idx = 1:4
                % Run Two-Way FE Model
                [pre_stats, post_stats, full_stats] = runTwoWayFEModel(month_offset_sorted, current_anomalies, radius_idx, 5);
                
                % Store in matrix
                stats_matrix(radius_idx, 1:6) = struct2array(pre_stats);
                stats_matrix(radius_idx, 7:12) = struct2array(post_stats);
                stats_matrix(radius_idx, 13:18) = struct2array(full_stats);
            end
            
            % Store in final results (Anomaly indices 5-8)
            final_results{anomaly_idx + 4}(ndvi_range_idx, :) = reshape(stats_matrix', 1, []);
        end
    end
    
    % Save 8 CSV files
    % Format depth and magnitude strings
    depth_str = num2str(round(EQ_depth * 10)); % Multiply depth by 10 and round
    mag_str = strrep(num2str(EQ_mag * 10, '%.0f'), '.', ''); % Multiply magnitude by 10 and remove decimal point
    land_prop_str = num2str(land_prop);
    % Modify output filenames
    output_names = {
        [num2str(EQ_date) '_' depth_str '_' mag_str '_1_' land_prop_str '_Standardized_Anomaly_Type1_Original_Monthly_Mean.csv'], ...
        [num2str(EQ_date) '_' depth_str '_' mag_str '_1_' land_prop_str '_Standardized_Anomaly_Type2_Below_Mean_Daily_Mean.csv'], ...
        [num2str(EQ_date) '_' depth_str '_' mag_str '_1_' land_prop_str '_Standardized_Anomaly_Type3_STL_Residual_Monthly_Mean.csv'], ...
        [num2str(EQ_date) '_' depth_str '_' mag_str '_1_' land_prop_str '_Standardized_Anomaly_Type4_STL_Negative_Residual_Mean.csv'], ...
        [num2str(EQ_date) '_' depth_str '_' mag_str '_1_' land_prop_str '_FE_Model_Type1_Original_Monthly_Mean.csv'], ...
        [num2str(EQ_date) '_' depth_str '_' mag_str '_1_' land_prop_str '_FE_Model_Type2_Below_Mean_Daily_Mean.csv'], ...
        [num2str(EQ_date) '_' depth_str '_' mag_str '_1_' land_prop_str '_FE_Model_Type3_STL_Residual_Monthly_Mean.csv'], ...
        [num2str(EQ_date) '_' depth_str '_' mag_str '_1_' land_prop_str '_FE_Model_Type4_STL_Negative_Residual_Mean.csv'] ...
    };
    
    % Create column names (modified for negative anomaly metrics)
    radius_names = {'R1_0_200km', 'R2_200_300km', 'R3_300_400km', 'R4_400_500km'};
    time_periods = {'Pre_EQ', 'Post_EQ', 'Full_Period'};
    stat_metrics = {'Mean_Negative_Anomaly', 'Cumulative_Negative_Intensity', ...
                   'Proportion_Significant_Months', 'Max_Negative_Anomaly', ...
                   'Mean_Anomaly', 'Max_Consecutive_Negative_Months'};
    
    % Generate complete column names
    column_names = {};
    for radius_idx = 1:4
        for time_idx = 1:3
            for stat_idx = 1:6
                col_name = sprintf('%s_%s_%s', radius_names{radius_idx}, ...
                                  time_periods{time_idx}, stat_metrics{stat_idx});
                column_names = [column_names, col_name];
            end
        end
    end
    
    % Save each CSV file
    for file_idx = 1:8
        % Create output table
        output_table = array2table(final_results{file_idx});
        output_table.Properties.VariableNames = column_names;
        
        % Add row names (NDVI Ranges) - Now 5 rows
        ndvi_ranges = {'NDVI_All';          % Row 1: No NDVI distinction
                      'NDVI_0_0.2';        % Row 2: Original 1st range
                      'NDVI_0.2_0.4';      % Row 3: Original 2nd range  
                      'NDVI_0.4_0.6';      % Row 4: Original 3rd range
                      'NDVI_0.6_1.0'};     % Row 5: Original 4th range
        output_table.NDVI_Range = ndvi_ranges;
        
        % Reorder columns, placing NDVI_Range first
        output_table = movevars(output_table, 'NDVI_Range', 'Before', 1);
        
        % Write table to CSV file at specified path
        output_filename = fullfile(output_folder, output_names{file_idx});
        writetable(output_table, output_filename);
        
    end
    
    fprintf('Saved all 8 CSV files for earthquake %d/%d successfully!\n', ij, height(eq_query));
    
    % Clear large matrices to free memory
    clear SSM_data_3D
end
