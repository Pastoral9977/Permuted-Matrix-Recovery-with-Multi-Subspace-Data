function [Inliers_id, Outliers_id, value,validation] ...
    = Outlier_Detect(X_tilde, detect_flag, alpha, validation, affine)

    [D,N] = size(X_tilde);
    if (nargin < 5)
        affine = false;
    end

    % "validation = true" means SoF scenario.
    if (nargin < 4) || (validation == -1)
        if (0.5*D < N)
            validation = true;
        else
            validation = false;
        end
    end

    if (nargin < 3) || (alpha == -1)
        alpha = 1e10;
    end

    if (nargin < 2)
        detect_flag = 0;
    end

    if (alpha <= 100)
        validation = false;
    end
    if (validation == false)
        alpha = computeLambda2_mat(X_tilde);
    end

    [Inliers_id, Outliers_id, value] = gank(X_tilde, affine, alpha, validation, detect_flag);

    invalue = value(Inliers_id);
    outvalue = value(Outliers_id);

    if (alpha >= 1e5)
        if (validation == true && prctile(invalue,75) >= prctile(outvalue,25))...
                || (validation == false && prctile(invalue, 25) <= prctile(outvalue,75))
            alpha = 10;
            validation = false;
            [Inliers_id, Outliers_id, ~] = gank(X_tilde, affine, alpha, validation, detect_flag);
        end
        elseif (alpha <= 100)
            if (prctile(invalue, 25) <= prctile(outvalue,75))
                alpha = 1e10;
                validation = true;
                [Inliers_id, Outliers_id, ~] = gank(X_tilde, affine, alpha, validation, detect_flag);
            end
    end    

    str = ['alpha = ' num2str(alpha) ', validation = '  num2str(validation)] ;
    disp(str);
end
 


function [Inliers_id, Outliers_id, value] = gank(X_tilde, affine, alpha, validation, detect_flag)

    if (detect_flag == 0)
        [Inliers_id, Outliers_id, value] = outlier_detect(X_tilde, affine, alpha, validation);

    elseif(detect_flag == 1)
        [Inliers_id, ~, value] = outlier_detect(X_tilde, affine, alpha, validation);
        disp(['Inliers_id = ' Inliers_id])
        [~,Outliers1_id] = outlier_detect(X_tilde(:,Inliers_id), affine, alpha, validation);
        Inliers_id = Inliers_id(setdiff(1:length(Inliers_id),Outliers1_id));
        Outliers_id = setdiff(1:size(X_tilde,2),Inliers_id)';

    elseif(detect_flag == 2)
        [Inliers_id, ~, value] = outlier_detect(X_tilde, affine, alpha, validation);
        disp(['Inliers_id = ' Inliers_id])
        [~,Outliers1_id] = outlier_detect(X_tilde(:,Inliers_id), affine, alpha, validation);
        Inliers_id = Inliers_id(setdiff(1:length(Inliers_id),Outliers1_id));
        Outliers_id = setdiff(1:size(X_tilde,2),Inliers_id);
        [Inliers2_id,~] = outlier_detect(X_tilde(:,Outliers_id), affine, alpha, validation);
        Outliers_id = Outliers_id(setdiff(1:length(Outliers_id),Inliers2_id));
        Inliers_id = setdiff(1:size(X_tilde,2),Outliers_id);
    end

end




function [Inliers_id, Outliers_id, value] = outlier_detect(X_tilde, affine, alpha, validation)
    Z = admmLasso_mat_func(X_tilde, affine, alpha);
    save CMat Z;
    value = vecnorm(Z,1);
    value = reshape(value, length(value), 1);
    if (validation)
        T = kmeans(value,2);
        I = find(T == 1); J = find(T == 2);
        I_mean = mean(value(I));J_mean = mean(value(J));
            if (I_mean > J_mean)
                Inliers_id = J;Outliers_id = I;
            else
                Inliers_id = I;Outliers_id = J;
            end
    else
        Outliers_id = find(value' < 1e-9);
        Inliers_id = setdiff(1:size(X_tilde,2),Outliers_id);
    end
end







