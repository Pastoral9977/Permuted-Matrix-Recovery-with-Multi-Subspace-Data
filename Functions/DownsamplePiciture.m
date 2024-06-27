% This function aims at downsampling a M*N picture to a m*n picture.
% Input: VECTORs of the original pictures
% Output: VECTORs of the downsampled pictures

function output = DownsamplePiciture(input, M, N, m, n)
    % Downsampling factors
    s = M / m;
    t = N / n;
    p = fix(s / 2); % �²���ÿ��ʱ����λ��ƫ������
    q = fix(t / 2); % �²���ʱȡ�е���λ��ƫ������

    % Initialize the output matrix
    output = zeros(m * n, size(input, 2));

    for k = 1:size(input, 2)
        % Reshape each column vector to MxN matrix
        X = reshape(input(:, k), M, N);

        % Initialize downsampled image matrix
        Y = zeros(m, n);

        for i = 1:n
            % Downsample each column of the image
            Y(:, i) = downsample(X(:, q + (i - 1) * t), s, p);
        end
        
        % Reshape downsampled image back to column vector
        output(:, k) = Y(:);
    end
end
