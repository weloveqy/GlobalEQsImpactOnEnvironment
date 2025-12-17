%% 不区分NDVI+区分DNVI  计算SSM数据；SSM负异常
clear
clc

% 读取文件
lst_filepath = "F:\GLDAS_SSM\"; 
ndvi_filepath = "F:\NDVI_mat\"; % NDVI filepath
output_folder = 'F:\Output_tables_figs';
eq_csvdata = "F:\EQs_anomaly_plot\EQs_query_screened_with_land_info.csv";
eq_query = readtable(eq_csvdata);

% 定义半径范围
radius_ranges = {
    '0-200km Circles', 0, 200;
    '200-300km Rings', 200, 300;
    '300-400km Rings', 300, 400;
    '400-500km Rings', 400, 500;
    '0-1000km Circles', 0, 1000;
    '0-1200km Circles', 0, 1200 
};

% 初始化网格和距离计算
resolution = 0.05; 
rows = 3600;
cols = 7200;
R = 6371; 
lat_grid = linspace(90 - resolution/2, -90 + resolution/2, rows);
lon_grid = linspace(-180 + resolution/2, 180 - resolution/2, cols);
[LON, LAT] = meshgrid(lon_grid, lat_grid);

% process loop 
for ij= 3:height(eq_query)
     
    % 参数设置与时间范围 
    EQ_time=eq_query.time(ij);
    EQ_date = strrep(EQ_time{1}(1:10), '-', '');
    EQ_year = str2num(EQ_date(1:4));
    EQ_month = str2num(EQ_date(5:6));
    EQ_lat = eq_query.latitude(ij); 
    EQ_lon = eq_query.longitude(ij);
    EQ_depth = eq_query.depth(ij);
    EQ_mag = eq_query.mag(ij);
    [dataset_start, dataset_end] = set_dataset_dates(EQ_date); % STL start and end input
    % 读取R1_Land_Prop值
    land_prop = eq_query.R1_Land_Prop(ij);

    % 检查条件：如果小于40，跳过当前循环；阈值为经验值，用于保留发生在洲际海岸带的地震，如南美智利
    if land_prop < 40 || EQ_lon > 157
        continue; % 跳入下一个循环
    end
    fprintf('==== 正在处理第%d/%d条地震数据，R1陆地占比: %d%% ====\n', ij, height(eq_query), land_prop);
    
    baseline_length = 48; % 48个月基线期, 24月影响期
    baseline_start = datetime(EQ_year, EQ_month, 1) - calmonths(baseline_length + 12);
    baseline_end = dateshift(datetime(EQ_year, EQ_month, 1) - calmonths(13), 'end', 'month');
    impact_start = datetime(EQ_year, EQ_month, 1) - calmonths(12);
    impact_end = dateshift(datetime(EQ_year, EQ_month, 1) + calmonths(11), 'end', 'month');
            
    dlat = deg2rad(LAT - EQ_lat);
    dlon = deg2rad(LON - EQ_lon);
    a = sin(dlat/2).^2 + cos(deg2rad(EQ_lat)) * cos(deg2rad(LAT)) .* sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distances = R * c; 
    
    % 创建1200km掩膜用于确定计算范围
    mask_1200km = distances <= 1200;
    
    % 找到包含1200km范围的最小边界框
    [row_indices, col_indices] = find(mask_1200km);
    if isempty(row_indices)
        fprintf('警告：地震位置 %f, %f 周围1200km范围内无有效数据，跳过处理\n', EQ_lat, EQ_lon);
        continue;
    end
    
    min_row = min(row_indices);
    max_row = max(row_indices);
    min_col = min(col_indices);
    max_col = max(col_indices);
    
    % 提取子区域的范围
    sub_rows = max_row - min_row + 1;
    sub_cols = max_col - min_col + 1;
    
    % 创建距离掩膜（使用子区域范围）
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
    
    % 获取所有需要处理的NDVI mat文件
    all_files = dir(fullfile(ndvi_filepath, 'NDVI_*.mat'));
    all_dates = [];
    
    for i = 1:length(all_files)
        filename = all_files(i).name;
        % 从文件名提取日期信息：NDVI_YYYY_DDD.mat
        date_parts = strsplit(filename(6:end-4), '_'); % 去掉"NDVI_"和".mat"
        if length(date_parts) == 2
            yeaR = str2double(date_parts{1});
            DOY = str2double(date_parts{2});
            date_val = datetime(yeaR, 1, 1) + days(DOY - 1);
            all_dates = [all_dates; date_val];
        end
    end
    
    % 计算NDVI掩膜
    % fprintf('开始计算NDVI时间掩膜...\n');
    ndvi_mask_start = datetime(EQ_year, EQ_month, 1) - calmonths(12);
    ndvi_mask_end = datetime(EQ_year, EQ_month, 1) + calmonths(12) - days(1);
    
    ndvi_mask_indices = (all_dates >= ndvi_mask_start) & (all_dates <= ndvi_mask_end);
    ndvi_mask_files = all_files(ndvi_mask_indices);
    ndvi_mask_dates = all_dates(ndvi_mask_indices);
    
    % 初始化NDVI累积数组（使用子区域范围）
    ndvi_accumulated = zeros(sub_rows, sub_cols);
    ndvi_count = zeros(sub_rows, sub_cols);
    fprintf('正在处理NDVI文件 %s ~ %s\n',ndvi_mask_start, ndvi_mask_end);

    for i = 1:length(ndvi_mask_files)
        filename = fullfile(ndvi_filepath, ndvi_mask_files(i).name);
        % fprintf('处理NDVI文件: %s\n', ndvi_mask_files(i).name);
        
        try
            % 加载.mat文件中的NDVI数据
            data = load(filename);
            NDVI_data = data.processed_data; 
            
            % 提取子区域
            NDVI_sub = NDVI_data(min_row:max_row, min_col:max_col);
                   
            spatial_mask = sub_distances <= 1000;
            valid_mask = spatial_mask & NDVI_sub ~= 0 & ~isnan(NDVI_sub) & NDVI_sub > -1;
            
            ndvi_accumulated(valid_mask) = ndvi_accumulated(valid_mask) + NDVI_sub(valid_mask);
            ndvi_count(valid_mask) = ndvi_count(valid_mask) + 1;
            
        catch ME
            fprintf('处理NDVI文件失败: %s, 错误: %s\n', filename, ME.message);
            continue;
        end
    end
    
    % NDVI范围设置 - 增加不区分NDVI的情况
    NDVI_ranges = [0, 1.0;  % 新增：不区分NDVI，覆盖全部范围
                   0, 0.2; 
                   0.2, 0.4; 
                   0.4, 0.6; 
                   0.6, 1.0];
    
    % 初始化结果存储 - 现在有5行（1行不区分NDVI + 4行区分NDVI）
    final_results = cell(8, 1); % 8种标准化异常指数
    for i = 1:8
        final_results{i} = zeros(5, 72); % 5行(1不区分NDVI + 4区分NDVI) × 72列(4半径×3区间×6指标)
    end
    
    % 提前读取所有SSM数据到一个三维矩阵中（使用子区域范围）
    fprintf('正在读取SSM数据 %s ~ %s\n', dataset_start, dataset_end);
    
    % 获取所有需要处理的SSM mat文件
    lst_all_files = dir(fullfile(lst_filepath, 'SSM_*.mat'));
    lst_all_dates = [];
    
    for i = 1:length(lst_all_files)
        filename = lst_all_files(i).name;
        % 从文件名提取日期信息：SSM_YYYY_DDD.mat
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
    
    % 初始化三维矩阵存储SSM数据（使用子区域范围）
    num_days = length(process_dates);
    SSM_data_3D = zeros(sub_rows, sub_cols, num_days);
    
    % 读取所有SSM数据（只读取子区域）
    for file_idx = 1:num_days
        filename = fullfile(lst_filepath, process_files(file_idx).name);
        
        try
            data = load(filename);
            % 只提取子区域的数据
            SSM_data_3D(:, :, file_idx) = data.SSM(min_row:max_row, min_col:max_col);
        catch ME
            fprintf('读取SSM文件失败: %s, 错误: %s\n', filename, ME.message);
            SSM_data_3D(:, :, file_idx) = NaN;
        end
    end
    
    % 循环处理5个NDVI范围（第1个是不区分NDVI的情况）
    for ndvi_range_idx = 1:5
        if ndvi_range_idx == 1
            fprintf('== 处理不区分NDVI范围 (1/5) ==\n');
        else
            fprintf('== 处理NDVI范围 %.1f-%.1f (%d/5) ==\n', ...
                NDVI_ranges(ndvi_range_idx, 1), NDVI_ranges(ndvi_range_idx, 2), ndvi_range_idx);
        end
        
        NDVI_threshold_low = NDVI_ranges(ndvi_range_idx, 1);
        NDVI_threshold_high = NDVI_ranges(ndvi_range_idx, 2);
    
        % 计算NDVI均值
        ndvi_mean = zeros(sub_rows, sub_cols);
        valid_pixels = ndvi_count > 0;
        ndvi_mean(valid_pixels) = ndvi_accumulated(valid_pixels) ./ ndvi_count(valid_pixels);
    
        % 创建NDVI阈值掩膜
        if ndvi_range_idx == 1
            % 不区分NDVI：使用所有有效像素
            ndvi_threshold_mask = valid_pixels;
        else
            % 区分NDVI：使用阈值范围
            ndvi_threshold_mask = (ndvi_mean >= NDVI_threshold_low) & (ndvi_mean <= NDVI_threshold_high) & valid_pixels;
        end
    
        % 初始化结果存储
        results = struct();
        results.daily_means = [];
        results.daily_stds = [];
        results.daily_counts = [];
        results.dates = [];
        
        % 处理每个SSM文件（直接从三维矩阵中获取数据）
        fprintf('正在处理SSM数据 %s ~ %s\n', dataset_start, dataset_end);

        for file_idx = 1:num_days
            current_date = process_dates(file_idx);
            
            % 直接从三维矩阵中获取SSM数据
            SSM_data = SSM_data_3D(:, :, file_idx);
            
            daily_means = zeros(1, 5);
            daily_stds = zeros(1, 5);
            daily_counts = zeros(1, 5);
            
            for radius_idx = 1:5 % 只处理前5个半径范围（0-200,200-300,300-400,400-500,0-1000）
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
    
        % 计算基线期月统计（原始SSM）
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
    
        % 计算影响期月统计（原始SSM）
        impact_mask = (results.dates >= impact_start) & (results.dates <= impact_end);
        impact_dates = results.dates(impact_mask);
        impact_means = results.daily_means(impact_mask, :);
        impact_years = year(impact_dates);
        impact_months = month(impact_dates);
        
        unique_month_pairs = unique([impact_years, impact_months], 'rows');
        
        % 初始化4种标准化异常指数（已修改为负异常）
        anomalies1 = []; % 第1种：原始月均值标准化异常
        anomalies2 = []; % 第2种：低于基线均值的日均值的月均值标准化异常
        anomalies3 = []; % 第3种：STL残差月均值标准化异常
        anomalies4 = []; % 第4种：STL负残差月均值标准化异常
        
        impact_month_offsets = [];
        
        for i = 1:size(unique_month_pairs, 1)
            year_val = unique_month_pairs(i, 1);
            month_val = unique_month_pairs(i, 2);
            
            month_mask = (impact_years == year_val) & (impact_months == month_val);
            
            monthly_mean1 = zeros(1, 5); % 原始月均值
            monthly_mean2 = zeros(1, 5); % 低于基线均值的日均值月均值
            monthly_anomaly1 = zeros(1, 5); % 第1种标准化异常
            monthly_anomaly2 = zeros(1, 5); % 第2种标准化异常
            
            for radius_idx = 1:5
                month_data = impact_means(month_mask, radius_idx);
                valid_data = month_data(~isnan(month_data));
                
                if ~isempty(valid_data)
                    % 原始月均值
                    monthly_mean1(radius_idx) = mean(valid_data);
                    
                    % 低于基线均值的日均值月均值
                    baseline_mean_val = baseline_monthly_means(month_val, radius_idx);
                    below_mean_data = valid_data(valid_data < baseline_mean_val);
                    
                    if ~isempty(below_mean_data)
                        monthly_mean2(radius_idx) = mean(below_mean_data);
                    else
                        monthly_mean2(radius_idx) = NaN;
                    end
                    
                    % 计算标准化异常
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
            
            % 计算月份偏移
            total_months = (year_val - EQ_year) * 12 + (month_val - EQ_month); % month offset
            month_offset = total_months;
            
            impact_month_offsets = [impact_month_offsets; month_offset];
            anomalies1 = [anomalies1; monthly_anomaly1];
            anomalies2 = [anomalies2; monthly_anomaly2];
        end
    
        % STL分解：对每个半径范围的逐日序列进行分解
        fprintf('开始执行STL分解...\n');
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
    
        % 计算基线期月统计（STL残差）
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
    
        % 计算影响期月统计（STL残差）
        impact_residuals = results.daily_residuals(impact_mask, :);
        
        % 计算第3和第4种标准化异常指数
        for i = 1:size(unique_month_pairs, 1)
            year_val = unique_month_pairs(i, 1);
            month_val = unique_month_pairs(i, 2);
            
            month_mask = (impact_years == year_val) & (impact_months == month_val);
            
            monthly_mean3 = zeros(1, 5); % STL残差月均值
            monthly_mean4 = zeros(1, 5); % STL负残差月均值
            monthly_anomaly3 = zeros(1, 5); % 第3种标准化异常
            monthly_anomaly4 = zeros(1, 5); % 第4种标准化异常
            
            for radius_idx = 1:5
                month_data = impact_residuals(month_mask, radius_idx);
                valid_data = month_data(~isnan(month_data));
                
                if ~isempty(valid_data)
                    % STL残差月均值
                    monthly_mean3(radius_idx) = mean(valid_data);
                    
                    % STL负残差月均值
                    baseline_residual_mean_val = baseline_residual_means(month_val, radius_idx);
                    below_residual_data = valid_data(valid_data < baseline_residual_mean_val);
                    
                    if ~isempty(below_residual_data)
                        monthly_mean4(radius_idx) = mean(below_residual_data);
                    else
                        monthly_mean4(radius_idx) = NaN;
                    end
                    
                    % 计算标准化异常
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
    
        % 按月份偏移排序
        [month_offset_sorted, sort_idx] = sort(impact_month_offsets);
        anomalies1_sorted = anomalies1(sort_idx, :);
        anomalies2_sorted = anomalies2(sort_idx, :);
        anomalies3_sorted = anomalies3(sort_idx, :);
        anomalies4_sorted = anomalies4(sort_idx, :);
    
        % 计算6个统计指标的函数（已修改为负异常统计）
        calculateStats = @(anomalies, mask) struct(...
            'MeanNegativeAnomaly', mean(anomalies(anomalies < 0 & mask)), ...
            'CumulativeNegativeIntensity', sum(anomalies(anomalies < 0 & mask)), ...
            'ProportionSignificant', sum(anomalies < 0 & mask) / sum(mask), ...
            'MaxNegativeAnomaly', min([anomalies(anomalies < 0 & mask); 0]), ...
            'MeanAnomaly', mean(anomalies(mask)), ...
            'MaxNegativeSequence', maxConsecutiveNegative(anomalies, mask) ...
        );
        
        % 为4种异常指数计算统计指标
        all_anomalies = {anomalies1_sorted, anomalies2_sorted, anomalies3_sorted, anomalies4_sorted};
        
        for anomaly_idx = 1:4
            current_anomalies = all_anomalies{anomaly_idx};
            stats_matrix = zeros(4, 18); % 4个半径 × 18个指标(3区间×6指标)
            
            for radius_idx = 1:4
                % 定义三个时间区间掩码
                pre_mask = (month_offset_sorted >= -11 & month_offset_sorted <= 0); % 震前12个月
                post_mask = (month_offset_sorted >= 0 & month_offset_sorted <= 11);  % 震后12个月
                full_mask = (month_offset_sorted >= -11 & month_offset_sorted <= 11); % 完整24个月
                
                % 计算三个时间区间的统计指标
                pre_stats = calculateStats(current_anomalies(:, radius_idx), pre_mask);
                post_stats = calculateStats(current_anomalies(:, radius_idx), post_mask);
                full_stats = calculateStats(current_anomalies(:, radius_idx), full_mask);
                
                % 存储到矩阵中
                stats_matrix(radius_idx, 1:6) = struct2array(pre_stats);
                stats_matrix(radius_idx, 7:12) = struct2array(post_stats);
                stats_matrix(radius_idx, 13:18) = struct2array(full_stats);
            end
            
            % 存储到最终结果中
            final_results{anomaly_idx}(ndvi_range_idx, :) = reshape(stats_matrix', 1, []);
        end
    
        % 双向固定效应模型分析（第5-8种异常指数）
        for anomaly_idx = 1:4
            current_anomalies = all_anomalies{anomaly_idx};
            stats_matrix = zeros(4, 18); % 4个半径 × 18个指标(3区间×6指标)
            
            for radius_idx = 1:4
                % 运行双向固定效应模型
                [pre_stats, post_stats, full_stats] = runTwoWayFEModel(month_offset_sorted, current_anomalies, radius_idx, 5);
                
                % 存储到矩阵中
                stats_matrix(radius_idx, 1:6) = struct2array(pre_stats);
                stats_matrix(radius_idx, 7:12) = struct2array(post_stats);
                stats_matrix(radius_idx, 13:18) = struct2array(full_stats);
            end
            
            % 存储到最终结果中（第5-8种异常指数）
            final_results{anomaly_idx + 4}(ndvi_range_idx, :) = reshape(stats_matrix', 1, []);
        end
    end
    
    % 保存8个CSV文件
    % 转换深度和震级格式
    depth_str = num2str(round(EQ_depth * 10)); % 深度乘以10后取整
    mag_str = strrep(num2str(EQ_mag * 10, '%.0f'), '.', ''); % 震级乘以10后去除小数点
    land_prop_str = num2str(land_prop);

    % 修改输出文件名
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
    
    % 创建列名（已修改为负异常指标）
    radius_names = {'R1_0_200km', 'R2_200_300km', 'R3_300_400km', 'R4_400_500km'};
    time_periods = {'Pre_EQ', 'Post_EQ', 'Full_Period'};
    stat_metrics = {'Mean_Negative_Anomaly', 'Cumulative_Negative_Intensity', ...
                   'Proportion_Significant_Months', 'Max_Negative_Anomaly', ...
                   'Mean_Anomaly', 'Max_Consecutive_Negative_Months'};
    
    % 生成完整的列名
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
    
    % 保存每个CSV文件
    for file_idx = 1:8
        % 创建输出表格
        output_table = array2table(final_results{file_idx});
        output_table.Properties.VariableNames = column_names;
        
        % 添加行名（NDVI范围）- 现在有5行
        ndvi_ranges = {'NDVI_All';          % 第1行：不区分NDVI
                      'NDVI_0_0.2';        % 第2行：原第1个范围
                      'NDVI_0.2_0.4';      % 第3行：原第2个范围  
                      'NDVI_0.4_0.6';      % 第4行：原第3个范围
                      'NDVI_0.6_1.0'};     % 第5行：原第4个范围
        output_table.NDVI_Range = ndvi_ranges;
        
        % 重新排列列，将NDVI_Range放在第一列
        output_table = movevars(output_table, 'NDVI_Range', 'Before', 1);
        
        % 写入CSV文件至指定路径
        output_filename = fullfile(output_folder, output_names{file_idx});
        writetable(output_table, output_filename);
        
    end
    
    fprintf('第%d/%d个地震的8个CSV文件均保存完成！\n', ij, height(eq_query));
    
    % 清除大型矩阵以释放内存
    clear SSM_data_3D
end

