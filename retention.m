function analyze_asrs_learning_retention_vr_vs_2d()
    % Load original and delayed CSV files, match participants,
    % analyze learning change and retention vs ASRS

    %% Select ORIGINAL files
    fprintf('=== Select ORIGINAL VR Group Files ===\n');
    [fn1, fp1] = uigetfile('*.csv', 'Select FIRST ORIGINAL VR CSV');
    if isequal(fn1,0), return; end
    [fn2, fp2] = uigetfile('*.csv', 'Select SECOND ORIGINAL VR CSV', fp1);
    if isequal(fn2,0), return; end
    
    fprintf('\n=== Select ORIGINAL 2D Group Files ===\n');
    [fn3, fp3] = uigetfile('*.csv', 'Select FIRST ORIGINAL 2D CSV', fp1);
    if isequal(fn3,0), return; end
    [fn4, fp4] = uigetfile('*.csv', 'Select SECOND ORIGINAL 2D CSV', fp3);
    if isequal(fn4,0), return; end
    
    %% Select DELAYED files
    fprintf('\n=== Select DELAYED VR Group Files ===\n');
    [dfn1, dfp1] = uigetfile('*.csv', 'Select FIRST DELAYED VR CSV', fp1);
    if isequal(dfn1,0), return; end
    [dfn2, dfp2] = uigetfile('*.csv', 'Select SECOND DELAYED VR CSV', dfp1);
    if isequal(dfn2,0), return; end
    
    fprintf('\n=== Select DELAYED 2D Group Files ===\n');
    [dfn3, dfp3] = uigetfile('*.csv', 'Select FIRST DELAYED 2D CSV', dfp1);
    if isequal(dfn3,0), return; end
    [dfn4, dfp4] = uigetfile('*.csv', 'Select SECOND DELAYED 2D CSV', dfp3);
    if isequal(dfn4,0), return; end
    
    %% Load files
    fprintf('\nLoading files...\n');
    T_vr_orig = combineFiles(fullfile(fp1,fn1), fullfile(fp2,fn2));
    T_2d_orig = combineFiles(fullfile(fp3,fn3), fullfile(fp4,fn4));
    T_vr_delay = combineDelayedFiles(fullfile(dfp1,dfn1), fullfile(dfp2,dfn2));
    T_2d_delay = combineDelayedFiles(fullfile(dfp3,dfn3), fullfile(dfp4,dfn4));
    
    %% Process groups
    fprintf('\n=== Processing VR Group ===\n');
    [vr_sum, vr_rl, vr_pl, vr_nl, vr_rr, vr_pr, vr_nr, vr_rd, vr_pd, vr_nd] = ...
        processGroupWithDelay(T_vr_orig, T_vr_delay, 'VR');
    
    fprintf('\n=== Processing 2D Group ===\n');
    [td_sum, td_rl, td_pl, td_nl, td_rr, td_pr, td_nr, td_rd, td_pd, td_nd] = ...
        processGroupWithDelay(T_2d_orig, T_2d_delay, '2D');
    
    %% Print comparisons
    fprintf('\n=== LEARNING CHANGE (Post-Pre) ===\n');
    fprintf('VR: N=%d, r=%.3f, p=%.4f\n', vr_nl, vr_rl, vr_pl);
    fprintf('2D: N=%d, r=%.3f, p=%.4f\n', td_nl, td_rl, td_pl);
    
    fprintf('\n=== RETENTION (Delayed-Post) ===\n');
    fprintf('VR: N=%d, r=%.3f, p=%.4f\n', vr_nr, vr_rr, vr_pr);
    fprintf('2D: N=%d, r=%.3f, p=%.4f\n', td_nr, td_rr, td_pr);
    
    fprintf('\n=== LEARNING OUTCOME 1 WEEK LATER (Delayed-Pre) ===\n');
    fprintf('VR: N=%d, r=%.3f, p=%.4f\n', vr_nd, vr_rd, vr_pd);
    fprintf('2D: N=%d, r=%.3f, p=%.4f\n', td_nd, td_rd, td_pd);
    
    %% Combine and analyze
    vr_sum.modality = repmat("VR", height(vr_sum), 1);
    td_sum.modality = repmat("2D", height(td_sum), 1);
    all_data = [vr_sum; td_sum];
    
    fprintf('\n=== HIGH vs LOW ASRS: LEARNING ===\n');
    analyzeLearningByASRS(all_data, 'learningChange', 'Learning Change');
    
    fprintf('\n=== HIGH vs LOW ASRS: RETENTION ===\n');
    analyzeLearningByASRS(all_data, 'retention', 'Retention');
    
    fprintf('\n=== HIGH vs LOW ASRS: LEARNING OUTCOME 1 WEEK LATER ===\n');
    analyzeLearningByASRS(all_data, 'overallDelayed', 'Learning Outcome 1 Week Later');
    
    %% Scatterplots
    makePlots(vr_sum, td_sum, vr_rl, vr_pl, td_rl, td_pl, vr_nl, td_nl, 'learningChange', 'Learning Change (Post-Pre)');
    makePlots(vr_sum, td_sum, vr_rr, vr_pr, td_rr, td_pr, vr_nr, td_nr, 'retention', 'Retention (Delayed-Post)');
    makePlots(vr_sum, td_sum, vr_rd, vr_pd, td_rd, td_pd, vr_nd, td_nd, 'overallDelayed', 'Learning Outcome 1 Week Later (Delayed-Pre)');
    
    %% Fisher tests
    fprintf('\n=== FISHER r-to-z TESTS ===\n');
    fprintf('Learning: '); runFisherTest(vr_rl, vr_nl, td_rl, td_nl);
    fprintf('Retention: '); runFisherTest(vr_rr, vr_nr, td_rr, td_nr);
    fprintf('Overall 1-week: '); runFisherTest(vr_rd, vr_nd, td_rd, td_nd);
    
    %% Save
    vr_sum.group = repmat("VR", height(vr_sum), 1);
    td_sum.group = repmat("2D", height(td_sum), 1);
    assignin('base', 'vr_summary', vr_sum);
    assignin('base', 'td_summary', td_sum);
    assignin('base', 'combined_summary', [vr_sum; td_sum]);
    fprintf('\nSaved: vr_summary, td_summary, combined_summary\n');
end

function T = combineFiles(f1, f2)
    T1 = readtable(f1, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    T2 = readtable(f2, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    cv = intersect(T1.Properties.VariableNames, T2.Properties.VariableNames, 'stable');
    T = [T1(:,cv); T2(:,cv)];
    fprintf('  Combined: %d + %d rows, %d columns\n', height(T1), height(T2), length(cv));
end

function T = combineDelayedFiles(f1, f2)
    % Read delayed files starting from row 4
    % Column R (18th column) = ID, Columns T-AC (20-29) = delayed test scores
    T1 = readtable(f1, 'TextType', 'string', 'VariableNamingRule', 'preserve', ...
                   'HeaderLines', 3, 'ReadVariableNames', false);
    T2 = readtable(f2, 'TextType', 'string', 'VariableNamingRule', 'preserve', ...
                   'HeaderLines', 3, 'ReadVariableNames', false);
    
    % Manually name the columns we need
    % Column 18 = R (ID), Columns 20-29 = T-AC (test scores)
    if width(T1) >= 29
        T1.Properties.VariableNames{18} = 'id';
        for i = 20:29
            T1.Properties.VariableNames{i} = sprintf('test_q%d', i-19);
        end
    end
    
    if width(T2) >= 29
        T2.Properties.VariableNames{18} = 'id';
        for i = 20:29
            T2.Properties.VariableNames{i} = sprintf('test_q%d', i-19);
        end
    end
    
    cv = intersect(T1.Properties.VariableNames, T2.Properties.VariableNames, 'stable');
    T = [T1(:,cv); T2(:,cv)];
    fprintf('  Combined delayed: %d + %d rows\n', height(T1), height(T2));
end

function [tbl, rl, pl, nl, rr, pr, nr, rd, pd, nd] = processGroupWithDelay(To, Td, gn)
    preV = "pre_test_q" + (1:10);
    postV = "post_test_q" + (1:10);
    asrsV = "asrs_q" + (1:18);
    
    c2m = @(t, v) str2double(string(table2array(t(:, v))));
    
    preM = c2m(To, preV);
    postM = c2m(To, postV);
    asrsM = c2m(To, asrsV);
    
    vr = any(~isnan(preM),2) | any(~isnan(postM),2) | any(~isnan(asrsM),2);
    preM = preM(vr,:); postM = postM(vr,:); asrsM = asrsM(vr,:);
    ids = string(To.id(vr));
    
    dv = string(Td.Properties.VariableNames);
    
    % Find ID column
    idColName = "";
    if ismember("id", dv)
        idColName = "id";
    elseif ismember("ID", dv)
        idColName = "ID";
    elseif ismember("Var18", dv)
        % Column R in Excel = column 18
        idColName = "Var18";
    else
        error('Could not find ID column in delayed data. Available columns: %s', strjoin(dv, ', '));
    end
    
    fprintf('Using ID column: %s\n', idColName);
    
    % Find delayed test columns - now they're named test_q1 through test_q10
    delV = "test_q" + (1:10);
    
    if all(ismember(delV, dv))
        fprintf('Using columns T-AC (test_q1-test_q10) for delayed test\n');
    else
        error('Could not find expected delayed test columns. Available: %s', strjoin(dv, ', '));
    end
    
    delM = c2m(Td, delV);
    idD = string(Td.(idColName));
    
    % Remove empty/NaN IDs and clean up whitespace
    idD = strtrim(idD);  % Remove leading/trailing whitespace
    validDelayedIDs = ~ismissing(idD) & strlength(idD) > 0 & strlength(idD) <= 4;  % Valid IDs are 3-4 digits
    idD = idD(validDelayedIDs);
    delM = delM(validDelayedIDs, :);
    
    fprintf('Found %d valid participants in delayed data (3-digit IDs)\n', length(idD));
    fprintf('Delayed IDs: %s\n', strjoin(idD, ', '));
    fprintf('Found %d participants in original data\n', length(ids));
    fprintf('Original IDs: %s\n', strjoin(ids, ', '));
    delS = nan(length(ids),1);
    matchCount = 0;
    noMatchIDs = {};
    for i = 1:length(ids)
        % Clean and trim both IDs for matching
        cleanOrigID = strtrim(ids(i));
        m = find(strcmp(idD, cleanOrigID));
        if ~isempty(m)
            delS(i) = sum(delM(m(1),:), 'omitnan');
            matchCount = matchCount + 1;
        else
            noMatchIDs{end+1} = char(cleanOrigID);
        end
    end
    
    fprintf('Successfully matched %d/%d participants with delayed data\n', matchCount, length(ids));
    if ~isempty(noMatchIDs)
        fprintf('WARNING: %d participants without delayed data:\n', length(noMatchIDs));
        fprintf('  Missing IDs: %s\n', strjoin(noMatchIDs, ', '));
    end
    
    pre = sum(preM, 2, 'omitnan');
    post = sum(postM, 2, 'omitnan');
    asrs = sum(asrsM, 2, 'omitnan');
    asrsA = sum(asrsM(:,1:6), 2, 'omitnan');
    asrsB = sum(asrsM(:,7:18), 2, 'omitnan');
    
    lc = post - pre;
    ret = delS - post;
    od = delS - pre;
    
    ag = repmat("Low", length(asrsA), 1);
    ag(asrsA >= 14) = "High";
    
    ac = repmat("", length(asrs), 1);
    ac(asrs <= 30) = "Low";
    ac(asrs >= 31 & asrs <= 39) = "Mild-Moderate";
    ac(asrs >= 40 & asrs <= 49) = "High";
    ac(asrs >= 50) = "Very High";
    
    [R1, P1] = corrcoef(asrs, lc, 'Rows', 'pairwise');
    rl = R1(1,2); pl = P1(1,2); nl = sum(~isnan(asrs) & ~isnan(lc));
    
    [R2, P2] = corrcoef(asrs, ret, 'Rows', 'pairwise');
    rr = R2(1,2); pr = P2(1,2); nr = sum(~isnan(asrs) & ~isnan(ret));
    
    [R3, P3] = corrcoef(asrs, od, 'Rows', 'pairwise');
    rd = R3(1,2); pd = P3(1,2); nd = sum(~isnan(asrs) & ~isnan(od));
    
    fprintf('N=%d, Learning r=%.3f p=%.4f, Retention r=%.3f p=%.4f\n', nl, rl, pl, rr, pr);
    
    tbl = table(ids, pre, post, delS, lc, ret, od, asrs, asrsA, asrsB, ag, ac, ...
        'VariableNames', {'participantID','preScore','postScore','delayedScore',...
        'learningChange','retention','overallDelayed','ASRS_total','ASRS_partA',...
        'ASRS_partB','ASRS_group','ASRS_total_category'});
end

function analyzeLearningByASRS(ad, m, ml)
    vl = ad(ad.modality=="VR" & ad.ASRS_group=="Low",:);
    vh = ad(ad.modality=="VR" & ad.ASRS_group=="High",:);
    tl = ad(ad.modality=="2D" & ad.ASRS_group=="Low",:);
    th = ad(ad.modality=="2D" & ad.ASRS_group=="High",:);
    
    mvl = mean(vl.(m),'omitnan'); mvh = mean(vh.(m),'omitnan');
    mtl = mean(tl.(m),'omitnan'); mth = mean(th.(m),'omitnan');
    
    svl = std(vl.(m),'omitnan')/sqrt(height(vl));
    svh = std(vh.(m),'omitnan')/sqrt(height(vh));
    stl = std(tl.(m),'omitnan')/sqrt(height(tl));
    sth = std(th.(m),'omitnan')/sqrt(height(th));
    
    fprintf('VR Low: M=%.2f SE=%.2f N=%d, VR High: M=%.2f SE=%.2f N=%d\n', mvl,svl,height(vl),mvh,svh,height(vh));
    fprintf('2D Low: M=%.2f SE=%.2f N=%d, 2D High: M=%.2f SE=%.2f N=%d\n', mtl,stl,height(tl),mth,sth,height(th));
    
    [~,pv,~,sv] = ttest2(vl.(m), vh.(m));
    [~,pt,~,st] = ttest2(tl.(m), th.(m));
    
    fprintf('VR: t(%d)=%.3f p=%.4f', sv.df, sv.tstat, pv);
    if pv<0.05, fprintf(' *\n'); else, fprintf('\n'); end
    fprintf('2D: t(%d)=%.3f p=%.4f', st.df, st.tstat, pt);
    if pt<0.05, fprintf(' *\n'); else, fprintf('\n'); end
    
    createBarChart(ml, mvl, mvh, mtl, mth, svl, svh, stl, sth, pv, pt, vl, vh, tl, th);
end

function makePlots(vs, ts, vr, vp, tr, tp, vn, tn, m, lbl)
    figure('Color','w','Position',[100,100,1200,500]);
    sgtitle(['ASRS vs ' lbl], 'FontSize', 14, 'FontWeight', 'bold');
    subplot(1,2,1);
    plotScatter(vs.ASRS_total, vs.(m), vr, vp, vn, 'VR', [0,0.45,0.74]);
    subplot(1,2,2);
    plotScatter(ts.ASRS_total, ts.(m), tr, tp, tn, '2D', [0.85,0.33,0.1]);
end

function plotScatter(x, y, r, p, n, ttl, clr)
    scatter(x, y, 60, clr, 'filled', 'MarkerFaceAlpha', 0.6); hold on;
    v = ~isnan(x) & ~isnan(y);
    if sum(v) >= 2
        pf = polyfit(x(v), y(v), 1);
        xf = linspace(min(x(v)), max(x(v)), 100);
        plot(xf, polyval(pf,xf), 'Color', clr, 'LineWidth', 2);
    end
    yline(0, '--k', 'LineWidth', 1); grid on;
    xlabel('ASRS Total Score', 'FontSize', 11);
    ylabel('Score Change', 'FontSize', 11);
    
    % Format p-value with significance stars
    pstr = sprintf('p=%.4f', p);
    if p < 0.001
        pstr = [pstr ' ***'];
    elseif p < 0.01
        pstr = [pstr ' **'];
    elseif p < 0.05
        pstr = [pstr ' *'];
    end
    
    title(sprintf('%s\nr=%.3f, %s\nN=%d', ttl, r, pstr, n), 'FontSize', 12);
end

function createBarChart(metricLabel, mvl, mvh, mtl, mth, svl, svh, stl, sth, pv, pt, vl, vh, tl, th)
    figure('Color','w','Position',[100,100,800,600]);
    
    means = [mvl, mvh; mtl, mth];
    errors = [svl, svh; stl, sth];
    
    b = bar(means, 'grouped');
    b(1).FaceColor = [0.3, 0.6, 0.9];
    b(2).FaceColor = [0.9, 0.4, 0.3];
    
    hold on;
    
    gw = 0.8/1.5;
    errorbar([1-gw/4,1+gw/4], [mvl,mvh], [svl,svh], 'k.', 'LineWidth',1.5);
    errorbar([2-gw/4,2+gw/4], [mtl,mth], [stl,sth], 'k.', 'LineWidth',1.5);
    
    % Add N labels above bars
    text(1-gw/4, mvl+svl+0.3, sprintf('N=%d', height(vl)), 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    text(1+gw/4, mvh+svh+0.3, sprintf('N=%d', height(vh)), 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    text(2-gw/4, mtl+stl+0.3, sprintf('N=%d', height(tl)), 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    text(2+gw/4, mth+sth+0.3, sprintf('N=%d', height(th)), 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    
    % Add significance markers
    ymax = max([mvl+svl, mvh+svh, mtl+stl, mth+sth]) + 1.5;
    
    % VR comparison
    if pv < 0.001
        text(1, ymax, '***', 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
    elseif pv < 0.01
        text(1, ymax, '**', 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
    elseif pv < 0.05
        text(1, ymax, '*', 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
    else
        text(1, ymax, 'n.s.', 'HorizontalAlignment','center','FontSize',11);
    end
    text(1, ymax+0.4, sprintf('p=%.4f', pv), 'HorizontalAlignment','center','FontSize',8);
    
    % 2D comparison
    if pt < 0.001
        text(2, ymax, '***', 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
    elseif pt < 0.01
        text(2, ymax, '**', 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
    elseif pt < 0.05
        text(2, ymax, '*', 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
    else
        text(2, ymax, 'n.s.', 'HorizontalAlignment','center','FontSize',11);
    end
    text(2, ymax+0.4, sprintf('p=%.4f', pt), 'HorizontalAlignment','center','FontSize',8);
    
    yline(0, '--k', 'LineWidth', 1);
    set(gca,'XTickLabel',{'VR','2D'});
    ylabel(metricLabel, 'FontSize', 12);
    xlabel('Modality', 'FontSize', 12);
    title([metricLabel ' by ASRS Group & Modality'], 'FontSize', 13, 'FontWeight', 'bold');
    legend({'Low ASRS (Part A < 14)','High ASRS (Part A ≥ 14)'},'Location','best','FontSize',10);
    grid on;
end

function runFisherTest(r1, n1, r2, n2)
    z1 = 0.5*log((1+r1)/(1-r1));
    z2 = 0.5*log((1+r2)/(1-r2));
    se = sqrt(1/(n1-3) + 1/(n2-3));
    z = (z1-z2)/se;
    p = 2*(1-normcdf(abs(z)));
    fprintf('VR r=%.3f vs 2D r=%.3f, z=%.3f, p=%.4f', r1, r2, z, p);
    if p<0.05, fprintf(' *\n'); else, fprintf('\n'); end
end