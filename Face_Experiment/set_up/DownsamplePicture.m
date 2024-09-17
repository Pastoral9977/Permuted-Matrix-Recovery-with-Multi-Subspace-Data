% This function aims at downsampling a M*N picture to a m*n picture.
% Input: VECTORs of the original pictures
% Output: VECTORs of the downsampled pictures

function output = DownsamplePicture(input, N, M, n, m)
    % Downsampling factors
    s = N / n;
    t = M / m;
    p = fix(s / 2); 
    q = ceil(t / 2); 

    % Initialize the output matrix
    output = zeros(n * m, size(input, 2));

    for k = 1:size(input, 2)
        % Reshape each column vector to MxN matrix
        X = reshape(input(:, k), N, M);

        % Initialize downsampled image matrix
        Y = zeros(n, m);

        for i = 1:m
            % Downsample each column of the image
            Y(:, i) = downsample(X(:, q + (i - 1) * t), s, p);
        end
        
        % Reshape downsampled image back to column vector
        output(:, k) = Y(:);
    end
end
