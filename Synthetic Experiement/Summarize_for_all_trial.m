

% 计算所有 trial 的均值、中位数和方差
mean_missrate_out = mean(cat(3, RESULT.missrate_out_alltrial{:}), 3);
median_missrate_out = median(cat(3, RESULT.missrate_out_alltrial{:}), 3);
var_missrate_out = var(cat(3, RESULT.missrate_out_alltrial{:}), 0, 3);

mean_missrate_out_gtbases = mean(cat(3, RESULT.missrate_out_gtbases_alltrial{:}), 3);
median_missrate_out_gtbases = median(cat(3, RESULT.missrate_out_gtbases_alltrial{:}), 3);
var_missrate_out_gtbases = var(cat(3, RESULT.missrate_out_gtbases_alltrial{:}), 0, 3);

mean_recover_error = mean(cat(3, RESULT.recover_error_alltrial{:}), 3);
median_recover_error = median(cat(3, RESULT.recover_error_alltrial{:}), 3);
var_recover_error = var(cat(3, RESULT.recover_error_alltrial{:}), 0, 3);

mean_recover_error_gtbases = mean(cat(3, RESULT.recover_error_gtbases_alltrial{:}), 3);
median_recover_error_gtbases = median(cat(3, RESULT.recover_error_gtbases_alltrial{:}), 3);
var_recover_error_gtbases = var(cat(3, RESULT.recover_error_gtbases_alltrial{:}), 0, 3);

mean_missrate_in = mean(cat(3, RESULT.missrate_in_alltrial{:}), 3);
median_missrate_in = median(cat(3, RESULT.missrate_in_alltrial{:}), 3);
var_missrate_in = var(cat(3, RESULT.missrate_in_alltrial{:}), 0, 3);

mean_undetected_outlier_ratio = mean(cat(3, RESULT.undetected_outlier_ratio_alltrial{:}), 3);
median_undetected_outlier_ratio = median(cat(3, RESULT.undetected_outlier_ratio_alltrial{:}), 3);
var_undetected_outlier_ratio = var(cat(3, RESULT.undetected_outlier_ratio_alltrial{:}), 0, 3);

% 将均值、中位数和方差存储在 RESULT 中
RESULT.mean_missrate_out = array2table(mean_missrate_out, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.median_missrate_out = array2table(median_missrate_out, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.var_missrate_out = array2table(var_missrate_out, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.mean_missrate_out_gtbases = array2table(mean_missrate_out_gtbases, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.median_missrate_out_gtbases = array2table(median_missrate_out_gtbases, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.var_missrate_out_gtbases = array2table(var_missrate_out_gtbases, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.mean_recover_error = array2table(mean_recover_error, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.median_recover_error = array2table(median_recover_error, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.var_recover_error = array2table(var_recover_error, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.mean_recover_error_gtbases = array2table(mean_recover_error_gtbases, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.median_recover_error_gtbases = array2table(median_recover_error_gtbases, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.var_recover_error_gtbases = array2table(var_recover_error_gtbases, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.mean_missrate_in = array2table(mean_missrate_in, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.median_missrate_in = array2table(median_missrate_in, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.var_missrate_in = array2table(var_missrate_in, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.mean_undetected_outlier_ratio = array2table(mean_undetected_outlier_ratio, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.median_undetected_outlier_ratio = array2table(median_undetected_outlier_ratio, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));

RESULT.var_undetected_outlier_ratio = array2table(var_undetected_outlier_ratio, ...
    'VariableNames', arrayfun(@(x) sprintf('ShuffledRatio_%.2f', x), shuffled_ratio_set, 'UniformOutput', false), ...
    'RowNames',      arrayfun(@(x) sprintf('RankRatio_%.2f', x),     rn_ratio_set,       'UniformOutput', false));
