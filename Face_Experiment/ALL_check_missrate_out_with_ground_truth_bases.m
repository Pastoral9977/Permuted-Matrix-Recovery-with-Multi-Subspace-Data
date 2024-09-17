close all
clear 

%% Set the paths
addpath(genpath(pwd));
addpath(genpath('..\Functions'));

%% Set global values
global Involve_reconstructed_bases DATA LABEL BASES RESULT INIT PIC

DATA = struct();
LABEL = struct();
BASES = struct();
RESULT = struct();
INIT = struct();
PIC = struct();


%% Set initial parameters

start_person_ids = [1, 6, 11, 16, 21, 26]; % Different start_person_id values
num_groups_range = [2, 3, 5,  8,  10, 12]; % num_groups from 3 to 10
% num_groups_range = [5,8,10,12]; % num_groups from 3 to 10

% start_person_ids = 1;
% num_groups_range = 10;

outliers_num_each = 16; % number of outliers out of 64 points
rrank = 8; 
missing_ratio = 0.0;
shuffled_ratio = 0.4;

INIT.outliers_num_each = outliers_num_each; % number of outliers out of 64 points
INIT.missing_ratio = missing_ratio;
INIT.shuffled_ratio = shuffled_ratio;
INIT.rrank = rrank; 
INIT.seeds = [];

Involve_reconstructed_bases = true; % Set this to true to include reconstructed bases

RESULT.missrate_out_gt_results              = zeros(length(num_groups_range), length(start_person_ids));
RESULT.missrate_out_reconstructed_results   = zeros(length(num_groups_range), length(start_person_ids));
RESULT.missrate_in_reconstructed_results    = zeros(length(num_groups_range), length(start_person_ids));
RESULT.undetected_outlier_num_ratio_results = zeros(length(num_groups_range), length(start_person_ids));
RESULT.recovery_err_recon_results           = zeros(length(num_groups_range), length(start_person_ids));
RESULT.recovery_err_gt_results              = zeros(length(num_groups_range), length(start_person_ids));

for ng_idx = 1:length(num_groups_range)
    num_groups = num_groups_range(ng_idx);
    INIT.num_groups = num_groups; 
    DATA.num_groups = num_groups;

    for sp_idx = 1:length(start_person_ids)
        seed = randi(1000000);
        INIT.seeds = [INIT.seeds, seed];
        rng(seed)
        start_person_id = start_person_ids(sp_idx);
        INIT.start_person_id = start_person_id; 
        
        fprintf([repmat('=', 1, 80),'\n'])
        fprintf(...
            "Progress %d/%d: num_groups = %d, start_person_id = %d, seed = %d\n",...
            (ng_idx-1)*length(start_person_ids)+sp_idx,...
            length(num_groups_range)*length(start_person_ids),...
            num_groups, ...
            start_person_id,...
            seed...
            )
        fprintf([repmat('=', 1, 80),'\n'])
        
        
        %% Prepare data and label
        run Data_Loading.m
        run Initialization.m
        run Data_Construction.m

        LABEL.allpoints_class_gt = zeros(N, 1);
        for ii = 1:DATA.num_groups
            LABEL.allpoints_class_gt(sum(nn(1:ii))+1:sum(nn(1:ii+1))) = ii;
        end

        %% Bases construction
        BASES.U_gt = ground_truth_bases_construction(DATA.X_gt, LABEL.allpoints_class_gt, INIT.rrank);

        if Involve_reconstructed_bases
            run Inliers_Detection.m
            run Basis_Reconstrution.m % obtain BASES.U_esti
            RESULT.missrate_in_reconstructed_results(ng_idx, sp_idx) = missrate_in;
        end

        %% Outlier classification   
        tic
        fprintf('Outlier classification...\n')

        % With ground truth bases
        fprintf('\twith ground truth bases...\n')
        [DATA.X_solved_gtbases, LABEL.outlier_class_solved_gtbases] = Solve_partial_LSR(DATA.X_tilde, LABEL.outliers_ID, BASES.U_gt);

        % With reconstructed bases
        if Involve_reconstructed_bases
            fprintf('\twith reconstructed bases...\n')
            [DATA.X_solved_reconstructedbases, LABEL.outlier_class_solved_reconstructedbases] = Solve_partial_LSR(DATA.X_tilde, LABEL.Outliers_id, BASES.U_esti);
        end
        fprintf('\tOver, took %0.2fs\n', toc)


       %% Assessment
        % For ground truth bases
        LABEL.allpoints_class_solved_gtbases = zeros(DATA.N, 1);
        for ii = 1:DATA.N
            if (ismember(ii, LABEL.outliers_ID) == 1)
                LABEL.allpoints_class_solved_gtbases(ii) = LABEL.outlier_class_solved_gtbases(LABEL.outliers_ID == ii);
            end
        end

        % For reconstructed bases
        LABEL.allpoints_class_solved_reconstructedbases = zeros(DATA.N, 1);
        if Involve_reconstructed_bases
            for ii = 1:DATA.N
                if ismember(ii, LABEL.Outliers_id) == 1
                    LABEL.allpoints_class_solved_reconstructedbases(ii) = LABEL.outlier_class_solved_reconstructedbases(LABEL.Outliers_id == ii);
                end
            end
        end
        
        % Calculate missrate for ground truth bases
        missrate_out_gt = sum(LABEL.allpoints_class_solved_gtbases(LABEL.outliers_ID) ~= LABEL.allpoints_class_gt(LABEL.outliers_ID)) / length(LABEL.outliers_ID);
        recovery_err_gt = Evaluate(DATA.X_solved_gtbases,DATA.X_gt);
        RESULT.missrate_out_gt_results(ng_idx, sp_idx) = missrate_out_gt;
        RESULT.recovery_err_gt_results(ng_idx, sp_idx) = recovery_err_gt;

        % Calculate missrate for reconstructed bases
        if Involve_reconstructed_bases
            undetected_outlier_num = length(setdiff(LABEL.outliers_ID, LABEL.Outliers_id));
            detected_true_Outliers_id = setdiff(LABEL.Outliers_id, LABEL.inliers_ID);
            missclassified_outlier_num = sum(LABEL.allpoints_class_solved_reconstructedbases(detected_true_Outliers_id) ~= LABEL.allpoints_class_gt(detected_true_Outliers_id));
            missrate_out_reconstructed = (missclassified_outlier_num) / (length(LABEL.outliers_ID)-undetected_outlier_num);
            RESULT.missrate_out_reconstructed_results(ng_idx, sp_idx) = missrate_out_reconstructed;
            undetected_outlier_num_ratio = undetected_outlier_num/length(LABEL.outliers_ID);
            RESULT.undetected_outlier_num_ratio_results(ng_idx, sp_idx) = undetected_outlier_num_ratio;
            recovery_err_recon = Evaluate(DATA.X_solved_reconstructedbases,DATA.X_gt);
            RESULT.recovery_err_recon_results(ng_idx, sp_idx) = recovery_err_recon;
        end
        fprintf('In this case, missrate_out_gt = %.3f, missrate_out_reconstructed = %.3f,\nwith undetected_outlier_num_ratio = %.3f\n\n',...
                missrate_out_gt, missrate_out_reconstructed, undetected_outlier_num_ratio)
    end
end

%% Calculate mean missrate
RESULT.mean_missrate_out_gt = mean(RESULT.missrate_out_gt_results, 2);
RESULT.mean_missrate_out_reconstructed = mean(RESULT.missrate_out_reconstructed_results, 2);
RESULT.mean_missrate_in_reconstructed = mean(RESULT.missrate_in_reconstructed_results, 2);
RESULT.mean_undetected_outlier_num_ratio_results = mean(RESULT.undetected_outlier_num_ratio_results, 2);

RESULT.median_missrate_out_gt = median(RESULT.missrate_out_gt_results, 2);
RESULT.median_missrate_out_reconstructed = median(RESULT.missrate_out_reconstructed_results, 2);
RESULT.median_missrate_in_reconstructed = median(RESULT.missrate_in_reconstructed_results, 2);
RESULT.median_undetected_outlier_num_ratio_results = median(RESULT.undetected_outlier_num_ratio_results, 2);

RESULT.mean_recovery_err_gt = mean(RESULT.recovery_err_gt_results, 2);
RESULT.mean_recovery_err_recon = mean(RESULT.recovery_err_recon_results, 2);
RESULT.median_recovery_err_gt = median(RESULT.recovery_err_gt_results, 2);
RESULT.median_recovery_err_recon = median(RESULT.recovery_err_recon_results, 2);

save('missrate_results.mat', 'RESULT');

%% Plot results
% Plot mean missrate out
figure;
plot(num_groups_range, RESULT.mean_missrate_out_gt, '-o', 'LineWidth', 2);
hold on;
plot(num_groups_range, RESULT.mean_missrate_out_reconstructed, '-x', 'LineWidth', 2);
hold off;
xlabel('Number of Groups');
ylabel('Mean Missrate Out');
title('Mean Missrate Out vs. Number of Groups');
legend('Ground Truth Bases', 'Reconstructed Bases');
grid on;
saveas(gcf, 'Mean_Missrate_Out_vs_Num_Groups.png'); % Save the mean missrate figure as a PNG file
savefig('Mean_Missrate_Out_vs_Num_Groups.fig'); % Save the mean missrate figure as a MATLAB figure file

% Plot median missrate out
figure;
plot(num_groups_range, RESULT.median_missrate_out_gt, '-o', 'LineWidth', 2);
hold on;
plot(num_groups_range, RESULT.median_missrate_out_reconstructed, '-x', 'LineWidth', 2);
hold off;
xlabel('Number of Groups');
ylabel('Median Missrate Out');
title('Median Missrate Out vs. Number of Groups');
legend('Ground Truth Bases', 'Reconstructed Bases');
grid on;
saveas(gcf, 'Median_Missrate_Out_vs_Num_Groups.png'); % Save the median missrate figure as a PNG file
savefig('Median_Missrate_Out_vs_Num_Groups.fig'); % Save the median missrate figure as a MATLAB figure file
