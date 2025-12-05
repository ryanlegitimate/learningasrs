function analyze_asrs_learning_vr_vs_2d()
    % ANALYZE_ASRS_LEARNING_VR_VS_2D
    % Load 2 CSV files for VR group and 2 for 2D group (counterbalanced),
    % combine them, compute pre/post scores and ASRS totals,
    % correlate ASRS with learning change for each modality,
    % and compare learning outcomes between VR and 2D.
    
    %% --- Select VR group files ---
    fprintf('=== Select VR Group Files ===\n');
    [fileName1, filePath1] = uigetfile({'*.csv','CSV files (*.csv)'}, ...
                                       'Select FIRST VR group CSV file');
    if isequal(fileName1,0)
        disp('File selection cancelled.');
        return;
    end
    
    [fileName2, filePath2] = uigetfile({'*.csv','CSV files (*.csv)'}, ...
                                       'Select SECOND VR group CSV file', ...
                                       filePath1);
    if isequal(fileName2,0)
        disp('File selection cancelled.');
        return;
    end
    
    %% --- Select 2D group files ---
    fprintf('\n=== Select 2D Group Files ===\n');
    [fileName3, filePath3] = uigetfile({'*.csv','CSV files (*.csv)'}, ...
                                       'Select FIRST 2D group CSV file', ...
                                       filePath1);
    if isequal(fileName3,0)
        disp('File selection cancelled.');
        return;
    end
    
    [fileName4, filePath4] = uigetfile({'*.csv','CSV files (*.csv)'}, ...
                                       'Select SECOND 2D group CSV file', ...
                                       filePath3);
    if isequal(fileName4,0)
        disp('File selection cancelled.');
        return;
    end
    
    %% --- Load and combine files ---
    vrFile1  = fullfile(filePath1, fileName1);
    vrFile2  = fullfile(filePath2, fileName2);
    tdFile1  = fullfile(filePath3, fileName3);
    tdFile2  = fullfile(filePath4, fileName4);
    
    fprintf('\nLoading files...\n');
    T_vr  = combineFiles(vrFile1, vrFile2);
    T_2d  = combineFiles(tdFile1, tdFile2);
    
    %% --- Process each group ---
    fprintf('\n=== Processing VR Group ===\n');
    [vr_summary, vr_r, vr_p, vr_n] = processGroup(T_vr, 'VR');
    
    fprintf('\n=== Processing 2D Group ===\n');
    [td_summary, td_r, td_p, td_n] = processGroup(T_2d, '2D');
    
    %% --- Compare groups ---
    fprintf('\n=== GROUP COMPARISON ===\n');
    fprintf('Total ASRS:\n');
    fprintf('  VR:  N=%d, r=%.3f, p=%.4f\n', vr_n, vr_r, vr_p);
    fprintf('  2D:  N=%d, r=%.3f, p=%.4f\n', td_n, td_r, td_p);
    
    %% --- High vs Low ASRS Analysis ---
    fprintf('\n=== HIGH vs LOW ASRS COMPARISON ===\n');
    
    % Add modality labels
    vr_summary.modality = repmat("VR", height(vr_summary), 1);
    td_summary.modality = repmat("2D", height(td_summary), 1);
    
    % Combine data
    all_data = [vr_summary; td_summary];
    
    % Calculate means for each combination of ASRS group × Modality
    vr_low = all_data(all_data.modality == "VR" & all_data.ASRS_group == "Low", :);
    vr_high = all_data(all_data.modality == "VR" & all_data.ASRS_group == "High", :);
    td_low = all_data(all_data.modality == "2D" & all_data.ASRS_group == "Low", :);
    td_high = all_data(all_data.modality == "2D" & all_data.ASRS_group == "High", :);
    
    mean_vr_low = mean(vr_low.learningChange, 'omitnan');
    mean_vr_high = mean(vr_high.learningChange, 'omitnan');
    mean_td_low = mean(td_low.learningChange, 'omitnan');
    mean_td_high = mean(td_high.learningChange, 'omitnan');
    
    se_vr_low = std(vr_low.learningChange, 'omitnan') / sqrt(height(vr_low));
    se_vr_high = std(vr_high.learningChange, 'omitnan') / sqrt(height(vr_high));
    se_td_low = std(td_low.learningChange, 'omitnan') / sqrt(height(td_low));
    se_td_high = std(td_high.learningChange, 'omitnan') / sqrt(height(td_high));
    
    fprintf('\nMean Learning Change by Group:\n');
    fprintf('  VR Low ASRS:   M=%.2f, SE=%.2f, N=%d\n', mean_vr_low, se_vr_low, height(vr_low));
    fprintf('  VR High ASRS:  M=%.2f, SE=%.2f, N=%d\n', mean_vr_high, se_vr_high, height(vr_high));
    fprintf('  2D Low ASRS:   M=%.2f, SE=%.2f, N=%d\n', mean_td_low, se_td_low, height(td_low));
    fprintf('  2D High ASRS:  M=%.2f, SE=%.2f, N=%d\n', mean_td_high, se_td_high, height(td_high));
    
    % T-tests within each modality
    [~, p_vr, ~, stats_vr] = ttest2(vr_low.learningChange, vr_high.learningChange);
    [~, p_2d, ~, stats_2d] = ttest2(td_low.learningChange, td_high.learningChange);
    
    fprintf('\nWithin-Modality Comparisons (Low vs High ASRS):\n');
    fprintf('  VR:  t(%d) = %.3f, p = %.4f', stats_vr.df, stats_vr.tstat, p_vr);
    if p_vr < 0.05
        fprintf(' *\n');
    else
        fprintf('\n');
    end
    fprintf('       SD Low = %.2f, SD High = %.2f\n', ...
            std(vr_low.learningChange, 'omitnan'), std(vr_high.learningChange, 'omitnan'));
    fprintf('       Cohen''s d = %.3f\n', ...
            (mean_vr_high - mean_vr_low) / sqrt((std(vr_low.learningChange, 'omitnan')^2 + std(vr_high.learningChange, 'omitnan')^2) / 2));
    
    fprintf('  2D:  t(%d) = %.3f, p = %.4f', stats_2d.df, stats_2d.tstat, p_2d);
    if p_2d < 0.05
        fprintf(' *\n');
    else
        fprintf('\n');
    end
    fprintf('       SD Low = %.2f, SD High = %.2f\n', ...
            std(td_low.learningChange, 'omitnan'), std(td_high.learningChange, 'omitnan'));
    fprintf('       Cohen''s d = %.3f\n', ...
            (mean_td_high - mean_td_low) / sqrt((std(td_low.learningChange, 'omitnan')^2 + std(td_high.learningChange, 'omitnan')^2) / 2));
    
    % 2-way ANOVA: Modality × ASRS Group
    fprintf('\n2-way ANOVA: Modality × ASRS Group\n');
    [p_anova, tbl_anova, stats_anova] = anovan(all_data.learningChange, ...
        {all_data.modality, all_data.ASRS_group}, ...
        'model', 'interaction', ...
        'varnames', {'Modality', 'ASRS_Group'}, ...
        'display', 'off');
    
    fprintf('Main effect of Modality: F(%d,%d)=%.2f, p=%.4f\n', ...
        tbl_anova{2,3}, tbl_anova{5,3}, tbl_anova{2,6}, p_anova(1));
    fprintf('Main effect of ASRS Group: F(%d,%d)=%.2f, p=%.4f\n', ...
        tbl_anova{3,3}, tbl_anova{5,3}, tbl_anova{3,6}, p_anova(2));
    fprintf('Interaction (Modality × ASRS): F(%d,%d)=%.2f, p=%.4f', ...
        tbl_anova{4,3}, tbl_anova{5,3}, tbl_anova{4,6}, p_anova(3));
    
    if p_anova(3) < 0.05
        fprintf(' ***\n');
        fprintf('Significant interaction! The effect of ASRS differs by modality.\n');
    else
        fprintf('\n');
    end
    
    % Mean learning change comparison
    vr_mean_learning = mean(vr_summary.learningChange, 'omitnan');
    td_mean_learning = mean(td_summary.learningChange, 'omitnan');
    fprintf('\nMean Learning Change:\n');
    fprintf('  VR: %.2f\n', vr_mean_learning);
    fprintf('  2D: %.2f\n', td_mean_learning);
    
    % Independent samples t-test
    [~, p_ttest, ~, stats] = ttest2(vr_summary.learningChange, ...
                                     td_summary.learningChange);
    fprintf('\nIndependent t-test (VR vs 2D learning change):\n');
    fprintf('  t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_ttest);
    
    %% --- Bar chart: High vs Low ASRS across modalities ---
    figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    
    means = [mean_vr_low, mean_vr_high; mean_td_low, mean_td_high];
    errors = [se_vr_low, se_vr_high; se_td_low, se_td_high];
    
    b = bar(means, 'grouped');
    b(1).FaceColor = [0.3, 0.6, 0.9];  % Low ASRS
    b(2).FaceColor = [0.9, 0.4, 0.3];  % High ASRS
    
    hold on;
    
    % Add error bars
    ngroups = size(means, 1);
    nbars = size(means, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, means(:,i), errors(:,i), 'k.', 'LineWidth', 1.5);
    end
    
    % Get sample sizes
    n_vr_low = height(vr_low);
    n_vr_high = height(vr_high);
    n_2d_low = height(td_low);
    n_2d_high = height(td_high);
    
    % Add N labels on bars
    x_vr_low = 1 - groupwidth/4;
    x_vr_high = 1 + groupwidth/4;
    x_2d_low = 2 - groupwidth/4;
    x_2d_high = 2 + groupwidth/4;
    
    text(x_vr_low, means(1,1) + errors(1,1) + 0.3, sprintf('N=%d', n_vr_low), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    text(x_vr_high, means(1,2) + errors(1,2) + 0.3, sprintf('N=%d', n_vr_high), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    text(x_2d_low, means(2,1) + errors(2,1) + 0.3, sprintf('N=%d', n_2d_low), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    text(x_2d_high, means(2,2) + errors(2,2) + 0.3, sprintf('N=%d', n_2d_high), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    
    % Add significance markers
    ymax = max(means(:) + errors(:)) + 1.5;
    
    % VR significance
    x_vr = 1;
    if p_vr < 0.001
        text(x_vr, ymax, '***', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    elseif p_vr < 0.01
        text(x_vr, ymax, '**', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    elseif p_vr < 0.05
        text(x_vr, ymax, '*', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    else
        text(x_vr, ymax, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    text(x_vr, ymax + 0.5, sprintf('p=%.3f', p_vr), 'HorizontalAlignment', 'center', 'FontSize', 8);
    
    % 2D significance
    x_2d = 2;
    if p_2d < 0.001
        text(x_2d, ymax, '***', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    elseif p_2d < 0.01
        text(x_2d, ymax, '**', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    elseif p_2d < 0.05
        text(x_2d, ymax, '*', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    else
        text(x_2d, ymax, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    text(x_2d, ymax + 0.5, sprintf('p=%.3f', p_2d), 'HorizontalAlignment', 'center', 'FontSize', 8);
    
    % Add zero line
    yline(0, '--k', 'LineWidth', 1);
    
    set(gca, 'XTickLabel', {'VR', '2D'});
    ylabel('Mean Learning Change (Post - Pre)');
    xlabel('Modality');
    title('Learning Change by Modality and ASRS Level');
    legend({'Low ASRS (Part A < 14)', 'High ASRS (Part A ≥ 14)'}, 'Location', 'best');
    grid on;
    
    hold off;
    
    %% --- Bar chart 2: Total Score categories across modalities ---
    figure('Color', 'w', 'Position', [150, 150, 1000, 600]);
    
    % Calculate means for each Total Score category × Modality
    categories = ["Low", "Mild-Moderate", "High", "Very High"];
    vr_means = zeros(1, 4);
    vr_ses = zeros(1, 4);
    vr_ns = zeros(1, 4);
    td_means = zeros(1, 4);
    td_ses = zeros(1, 4);
    td_ns = zeros(1, 4);
    
    for i = 1:4
        cat = categories(i);
        
        % VR group
        vr_cat = all_data(all_data.modality == "VR" & all_data.ASRS_total_category == cat, :);
        if height(vr_cat) > 0
            vr_means(i) = mean(vr_cat.learningChange, 'omitnan');
            vr_ses(i) = std(vr_cat.learningChange, 'omitnan') / sqrt(height(vr_cat));
            vr_ns(i) = height(vr_cat);
        else
            vr_means(i) = NaN;
            vr_ses(i) = 0;
            vr_ns(i) = 0;
        end
        
        % 2D group
        td_cat = all_data(all_data.modality == "2D" & all_data.ASRS_total_category == cat, :);
        if height(td_cat) > 0
            td_means(i) = mean(td_cat.learningChange, 'omitnan');
            td_ses(i) = std(td_cat.learningChange, 'omitnan') / sqrt(height(td_cat));
            td_ns(i) = height(td_cat);
        else
            td_means(i) = NaN;
            td_ses(i) = 0;
            td_ns(i) = 0;
        end
    end
    
    % Create grouped bar chart
    means_total = [vr_means; td_means]';
    
    b = bar(means_total, 'grouped');
    b(1).FaceColor = [0, 0.4470, 0.7410];  % VR
    b(2).FaceColor = [0.8500, 0.3250, 0.0980];  % 2D
    
    hold on;
    
    % Add error bars
    ngroups = size(means_total, 1);
    nbars = size(means_total, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        if i == 1
            errorbar(x, vr_means, vr_ses, 'k.', 'LineWidth', 1.5);
        else
            errorbar(x, td_means, td_ses, 'k.', 'LineWidth', 1.5);
        end
    end
    
    % Add N labels
    ymax_total = max(max([vr_means; td_means] + [vr_ses; td_ses], [], 'omitnan')) + 0.5;
    for i = 1:ngroups
        x_vr = i - groupwidth/4;
        x_2d = i + groupwidth/4;
        
        if vr_ns(i) > 0
            text(x_vr, vr_means(i) + vr_ses(i) + 0.15, sprintf('N=%d', vr_ns(i)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
        
        if td_ns(i) > 0
            text(x_2d, td_means(i) + td_ses(i) + 0.15, sprintf('N=%d', td_ns(i)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    
    % Add zero line
    yline(0, '--k', 'LineWidth', 1);
    
    set(gca, 'XTickLabel', categories);
    ylabel('Mean Learning Change (Post - Pre)');
    xlabel('ASRS Total Score Category');
    title('Learning Change by ASRS Total Score Category and Modality');
    legend({'VR', '2D'}, 'Location', 'best');
    grid on;
    xtickangle(15);
    
    hold off;
    
    fprintf('\n=== TOTAL SCORE CATEGORIES BREAKDOWN ===\n');
    for i = 1:4
        fprintf('%s:\n', categories(i));
        fprintf('  VR:  M=%.2f, SE=%.2f, N=%d\n', vr_means(i), vr_ses(i), vr_ns(i));
        fprintf('  2D:  M=%.2f, SE=%.2f, N=%d\n', td_means(i), td_ses(i), td_ns(i));
    end
    
    %% --- Bar chart 3: Low ASRS Score group only ---
    figure('Color', 'w', 'Position', [200, 200, 600, 600]);
    
    % Get Low ASRS data (Total Score ≤ 30)
    vr_low_asrs = all_data(all_data.modality == "VR" & all_data.ASRS_total_category == "Low", :);
    td_low_asrs = all_data(all_data.modality == "2D" & all_data.ASRS_total_category == "Low", :);
    
    mean_vr_low_asrs = mean(vr_low_asrs.learningChange, 'omitnan');
    mean_td_low_asrs = mean(td_low_asrs.learningChange, 'omitnan');
    se_vr_low_asrs = std(vr_low_asrs.learningChange, 'omitnan') / sqrt(height(vr_low_asrs));
    se_td_low_asrs = std(td_low_asrs.learningChange, 'omitnan') / sqrt(height(td_low_asrs));
    n_vr_low_asrs = height(vr_low_asrs);
    n_td_low_asrs = height(td_low_asrs);
    
    % T-test comparing VR vs 2D in Low ASRS group
    [~, p_low_asrs, ~, stats_low_asrs] = ttest2(vr_low_asrs.learningChange, td_low_asrs.learningChange);
    
    % Create bar chart
    means_low = [mean_vr_low_asrs, mean_td_low_asrs];
    errors_low = [se_vr_low_asrs, se_td_low_asrs];
    
    b = bar(means_low);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0, 0.4470, 0.7410];  % VR
    b.CData(2,:) = [0.8500, 0.3250, 0.0980];  % 2D
    
    hold on;
    
    % Add error bars
    errorbar(1:2, means_low, errors_low, 'k.', 'LineWidth', 2, 'MarkerSize', 1);
    
    % Add N labels
    text(1, mean_vr_low_asrs + se_vr_low_asrs + 0.2, sprintf('N=%d', n_vr_low_asrs), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    text(2, mean_td_low_asrs + se_td_low_asrs + 0.2, sprintf('N=%d', n_td_low_asrs), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    
    % Add p-value
    ymax_low = max(means_low + errors_low) + 0.8;
    if p_low_asrs < 0.001
        text(1.5, ymax_low, '***', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    elseif p_low_asrs < 0.01
        text(1.5, ymax_low, '**', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    elseif p_low_asrs < 0.05
        text(1.5, ymax_low, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    else
        text(1.5, ymax_low, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    text(1.5, ymax_low + 0.3, sprintf('p=%.3f', p_low_asrs), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
    
    % Add zero line
    yline(0, '--k', 'LineWidth', 1.5);
    
    set(gca, 'XTickLabel', {'VR', '2D'});
    ylabel('Mean Learning Change (Post - Pre)', 'FontSize', 12);
    xlabel('Modality', 'FontSize', 12);
    title('Low ASRS Score', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);
    
    hold off;
    
    fprintf('\n=== LOW ASRS SCORE GROUP (Total ≤ 30) ===\n');
    fprintf('VR:  M=%.2f, SE=%.2f, N=%d\n', mean_vr_low_asrs, se_vr_low_asrs, n_vr_low_asrs);
    fprintf('2D:  M=%.2f, SE=%.2f, N=%d\n', mean_td_low_asrs, se_td_low_asrs, n_td_low_asrs);
    fprintf('t(%d) = %.3f, p = %.4f', stats_low_asrs.df, stats_low_asrs.tstat, p_low_asrs);
    if p_low_asrs < 0.05
        fprintf(' *\n');
    else
        fprintf('\n');
    end
    
    %% --- Individual group scatterplots ---
    figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
    
    % VR subplot
    subplot(1,2,1);
    plotGroup(vr_summary.ASRS_total, vr_summary.learningChange, vr_r, vr_p, 'VR', [0, 0.4470, 0.7410]);
    
    % 2D subplot
    subplot(1,2,2);
    plotGroup(td_summary.ASRS_total, td_summary.learningChange, td_r, td_p, '2D', [0.8500, 0.3250, 0.0980]);
    
    %% --- Combined correlation plot (all participants) ---
    figure('Color', 'w');
    
    % Plot VR group
    scatter(vr_summary.ASRS_total, vr_summary.learningChange, 60, ...
            [0, 0.4470, 0.7410], 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Plot 2D group
    scatter(td_summary.ASRS_total, td_summary.learningChange, 60, ...
            [0.8500, 0.3250, 0.0980], 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Overall correlation (both groups combined)
    all_asrs = [vr_summary.ASRS_total; td_summary.ASRS_total];
    all_learning = [vr_summary.learningChange; td_summary.learningChange];
    
    [R_all, P_all] = corrcoef(all_asrs, all_learning, 'Rows', 'pairwise');
    r_all = R_all(1,2);
    p_all = P_all(1,2);
    n_all = sum(~isnan(all_asrs) & ~isnan(all_learning));
    
    % Overall regression line
    mask_all = ~isnan(all_asrs) & ~isnan(all_learning);
    x_all = all_asrs(mask_all);
    y_all = all_learning(mask_all);
    if numel(x_all) >= 2
        p_fit = polyfit(x_all, y_all, 1);
        xfit = linspace(min(x_all), max(x_all), 100);
        yfit = polyval(p_fit, xfit);
        plot(xfit, yfit, 'k-', 'LineWidth', 2);
    end
    
    % Zero-change line
    yline(0, '--', 'Zero change', ...
          'LabelHorizontalAlignment', 'left', ...
          'LabelVerticalAlignment',   'bottom');
    
    grid on;
    xlabel('ASRS Total Score');
    ylabel('Learning Change (Post - Pre)');
    title(sprintf('Combined: ASRS vs Learning Change\n(r = %.2f, p = %.3f, N = %d)', ...
                  r_all, p_all, n_all));
    legend({'VR', '2D', 'Overall Fit'}, 'Location', 'best');
    
    % Adjust axes
    xlim([min(all_asrs)-1, max(all_asrs)+1]);
    yRange = range(all_learning);
    yMargin = max(1, 0.1 * yRange);
    ylim([min(all_learning)-yMargin, max(all_learning)+yMargin]);
    
    hold off;
    
    fprintf('\nOverall Correlation (both groups combined):\n');
    fprintf('  r = %.3f, p = %.4f, N = %d\n', r_all, p_all, n_all);
    
    %% --- Fisher's r-to-z test (compare correlations) ---
    fprintf('\n=== COMPARING CORRELATIONS (Fisher''s r-to-z) ===\n');
    
    % Total ASRS
    fprintf('\nTotal ASRS vs Learning:\n');
    z_vr = 0.5 * log((1 + vr_r) / (1 - vr_r));
    z_2d = 0.5 * log((1 + td_r) / (1 - td_r));
    se_diff = sqrt((1 / (vr_n - 3)) + (1 / (td_n - 3)));
    z_stat = (z_vr - z_2d) / se_diff;
    p_fisher = 2 * (1 - normcdf(abs(z_stat)));
    
    fprintf('  VR: r = %.3f, 2D: r = %.3f\n', vr_r, td_r);
    fprintf('  Fisher''s z = %.3f, p = %.4f', z_stat, p_fisher);
    if p_fisher < 0.05
        fprintf(' ***\n');
    else
        fprintf('\n');
    end
    
    %% --- Moderated regression (interaction test) ---
    fprintf('\n=== MODERATED REGRESSION (Interaction Test) ===\n');
    
    % Create modality dummy variable (0 = 2D, 1 = VR)
    modality_vr = ones(height(vr_summary), 1);
    modality_2d = zeros(height(td_summary), 1);
    modality_all = [modality_vr; modality_2d];
    
    % Prepare data
    all_asrs = [vr_summary.ASRS_total; td_summary.ASRS_total];
    all_learning = [vr_summary.learningChange; td_summary.learningChange];
    
    % Total ASRS interaction
    fprintf('\nTotal ASRS × Modality:\n');
    runModeratedRegression(all_asrs, modality_all, all_learning, 'Total ASRS');
    
    %% --- Save to workspace ---
    vr_summary.group = repmat("VR", height(vr_summary), 1);
    td_summary.group = repmat("2D", height(td_summary), 1);
    
    combined_summary = [vr_summary; td_summary];
    
    assignin('base', 'vr_summary', vr_summary);
    assignin('base', 'td_summary', td_summary);
    assignin('base', 'combined_summary', combined_summary);
    
    fprintf('\nSummary tables saved to workspace:\n');
    fprintf('  - vr_summary\n');
    fprintf('  - td_summary\n');
    fprintf('  - combined_summary\n');
end

%% --- HELPER FUNCTIONS ---

function T_combined = combineFiles(file1, file2)
    % Load both files and combine them
    T1 = readtable(file1, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    T2 = readtable(file2, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    
    % Get variable names from both tables
    vars1 = T1.Properties.VariableNames;
    vars2 = T2.Properties.VariableNames;
    
    % Find common variables
    commonVars = intersect(vars1, vars2, 'stable');
    
    if isempty(commonVars)
        error('No common column names found between the two files.');
    end
    
    % Keep only common columns from both tables
    T1_subset = T1(:, commonVars);
    T2_subset = T2(:, commonVars);
    
    % Combine vertically
    T_combined = [T1_subset; T2_subset];
    
    fprintf('  Combined %d rows from file 1 and %d rows from file 2\n', ...
            height(T1), height(T2));
    fprintf('  Using %d common columns\n', length(commonVars));
end

function [summaryTbl, r, pVal, nObs] = processGroup(T, groupName)
    % Process a single group's data
    
    %% --- Expected columns ---
    preVars   = "pre_test_q"  + (1:10);
    postVars  = "post_test_q" + (1:10);
    asrsVars  = "asrs_q"      + (1:18);
    asrsPartA = "asrs_q"      + (1:6);    % Part A: most predictive items
    asrsPartB = "asrs_q"      + (7:18);   % Part B: additional symptoms
    idVarName = "id";
    
    % Check presence of required columns
    requiredCols = [preVars, postVars, asrsVars, idVarName];
    missing = setdiff(requiredCols, string(T.Properties.VariableNames));
    if ~isempty(missing)
        error('Missing expected columns in %s group: %s', ...
              groupName, strjoin(missing, ', '));
    end
    
    %% --- Helper: convert table columns to numeric matrix ---
    colToMatrix = @(tbl, varNames) ...
        str2double(string(table2array(tbl(:, varNames))));
    
    %% --- Convert all rows ---
    preMatAll  = colToMatrix(T, preVars);
    postMatAll = colToMatrix(T, postVars);
    asrsMatAll = colToMatrix(T, asrsVars);
    
    %% --- Detect valid participant rows ---
    validRows = any(~isnan(preMatAll), 2) | ...
                any(~isnan(postMatAll), 2) | ...
                any(~isnan(asrsMatAll), 2);
    
    if nnz(validRows) == 0
        error('No valid participant rows in %s group.', groupName);
    end
    
    % Subselect matrices and ID column
    preMat   = preMatAll(validRows, :);
    postMat  = postMatAll(validRows, :);
    asrsMat  = asrsMatAll(validRows, :);
    idRaw    = T.(idVarName);
    idSub    = idRaw(validRows);
    
    %% --- Compute scores ---
    preScore       = sum(preMat,  2, 'omitnan');
    postScore      = sum(postMat, 2, 'omitnan');
    asrsTotal      = sum(asrsMat, 2, 'omitnan');
    asrsPartAScore = sum(asrsMat(:, 1:6), 2, 'omitnan');    % Part A (items 1-6)
    asrsPartBScore = sum(asrsMat(:, 7:18), 2, 'omitnan');   % Part B (items 7-18)
    learningChange = postScore - preScore;
    
    % Classify as Low (<14) or High (≥14) based on Part A
    asrsGroup = repmat("Low", length(asrsPartAScore), 1);
    asrsGroup(asrsPartAScore >= 14) = "High";
    
    % Classify by Total Score categories
    asrsTotalCategory = repmat("", length(asrsTotal), 1);
    asrsTotalCategory(asrsTotal <= 30) = "Low";
    asrsTotalCategory(asrsTotal >= 31 & asrsTotal <= 39) = "Mild-Moderate";
    asrsTotalCategory(asrsTotal >= 40 & asrsTotal <= 49) = "High";
    asrsTotalCategory(asrsTotal >= 50) = "Very High";
    
    %% --- Pearson correlation ---
    [R, P] = corrcoef(asrsTotal, learningChange, 'Rows', 'pairwise');
    r     = R(1,2);
    pVal  = P(1,2);
    nObs  = sum(~isnan(asrsTotal) & ~isnan(learningChange));
    
    fprintf('N = %d\n', nObs);
    fprintf('Pearson r (ASRS Total vs learning) = %.3f, p = %.4f\n', r, pVal);
    
    % Count participants in each ASRS group
    nLow = sum(asrsGroup == "Low");
    nHigh = sum(asrsGroup == "High");
    fprintf('Part A - Low ASRS (< 14): N = %d\n', nLow);
    fprintf('Part A - High ASRS (≥ 14): N = %d\n', nHigh);
    
    % Count participants in each Total Score category
    nTotalLow = sum(asrsTotalCategory == "Low");
    nTotalMild = sum(asrsTotalCategory == "Mild-Moderate");
    nTotalHigh = sum(asrsTotalCategory == "High");
    nTotalVeryHigh = sum(asrsTotalCategory == "Very High");
    fprintf('Total Score - Low (≤30): N = %d\n', nTotalLow);
    fprintf('Total Score - Mild-Moderate (31-39): N = %d\n', nTotalMild);
    fprintf('Total Score - High (40-49): N = %d\n', nTotalHigh);
    fprintf('Total Score - Very High (≥50): N = %d\n', nTotalVeryHigh);
    
    %% --- Create summary table ---
    summaryTbl = table(idSub, preScore, postScore, learningChange, asrsTotal, ...
                       asrsPartAScore, asrsPartBScore, asrsGroup, asrsTotalCategory, ...
        'VariableNames', {'participantID', 'preScore', 'postScore', ...
                          'learningChange', 'ASRS_total', 'ASRS_partA', ...
                          'ASRS_partB', 'ASRS_group', 'ASRS_total_category'});
end

function plotGroup(asrs_scores, learningChange, r, pVal, titleText, color)
    % Create scatterplot for a group
    
    nObs = sum(~isnan(asrs_scores) & ~isnan(learningChange));
    
    scatter(asrs_scores, learningChange, 60, color, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Regression line
    mask = ~isnan(asrs_scores) & ~isnan(learningChange);
    x = asrs_scores(mask);
    y = learningChange(mask);
    if numel(x) >= 2
        p = polyfit(x, y, 1);
        xfit = linspace(min(x), max(x), 100);
        yfit = polyval(p, xfit);
        plot(xfit, yfit, 'Color', color, 'LineWidth', 2);
    end
    
    % Zero-change line
    yline(0, '--k', 'Zero change', ...
          'LabelHorizontalAlignment', 'left', ...
          'LabelVerticalAlignment',   'bottom', ...
          'FontSize', 8);
    
    grid on;
    xlabel('ASRS Total Score');
    ylabel('Learning Change (Post - Pre)');
    title(sprintf('%s: ASRS vs Learning\n(r = %.2f, p = %.3f, N = %d)', ...
                  titleText, r, pVal, nObs));
    
    % Adjust axes
    if ~isempty(x)
        xlim([min(asrs_scores)-1, max(asrs_scores)+1]);
        yRange  = range(learningChange);
        yMargin = max(1, 0.1 * yRange);
        ylim([min(learningChange)-yMargin, max(learningChange)+yMargin]);
    end
    
    hold off;
end

function runModeratedRegression(predictor, modality, outcome, predictorName)
    % Run moderated regression and test interaction
    
    % Remove NaN cases
    valid_idx = ~isnan(predictor) & ~isnan(outcome) & ~isnan(modality);
    pred_clean = predictor(valid_idx);
    outcome_clean = outcome(valid_idx);
    mod_clean = modality(valid_idx);
    
    % Create interaction term
    interaction = pred_clean .* mod_clean;
    
    % Full model: outcome ~ predictor + modality + interaction
    X = [ones(length(pred_clean), 1), pred_clean, mod_clean, interaction];
    [b, ~, ~, ~, stats] = regress(outcome_clean, X);
    
    fprintf('   %s coefficient: b = %.3f\n', predictorName, b(2));
    fprintf('   Modality coefficient: b = %.3f\n', b(3));
    fprintf('   Interaction coefficient: b = %.3f\n', b(4));
    
    % Reduced model (no interaction)
    X_reduced = [ones(length(pred_clean), 1), pred_clean, mod_clean];
    [~, ~, ~, ~, stats_reduced] = regress(outcome_clean, X_reduced);
    
    % F-test for interaction
    SSR_reduced = stats_reduced(4);
    SSR_full = stats(4);
    df_diff = 1;
    df_error = length(pred_clean) - 4;
    
    F_interaction = ((SSR_reduced - SSR_full) / df_diff) / (SSR_full / df_error);
    p_interaction = 1 - fcdf(F_interaction, df_diff, df_error);
    
    fprintf('   Interaction F(%d,%d) = %.3f, p = %.4f', df_diff, df_error, F_interaction, p_interaction);
    if p_interaction < 0.05
        fprintf(' ***\n');
    else
        fprintf('\n');
    end
end