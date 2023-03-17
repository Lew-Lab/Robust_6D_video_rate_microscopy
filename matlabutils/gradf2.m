function g = gradf2(gamma, SMLM_img, b, FPSF, recovStruct)
%gradf computes the gradient of the negative Poisson log-likelihood
%w.r.t. molecular parameters at the current estimate gamma
%->-----
% input
%->-----
% gamma:        array (N*3,n_f) -current estimate of the molecular parameters
% SMLM_img:     array (M,n_f)   -input camera images
% b:            array (M,n_f)   -background

% N:                            -number of grid points
% M:                            -image size in number of camera pixels
% n_f:                          -number of frames
% upsamplef:                    -object space upsample factor
%----->-
% output
%----->-
% g:                            -array (N*3,n_f) gradient of the negative Poisson
%                                log-likelihood w.r.t. gamma

%% 1- global parameters

N = recovStruct.n_grid_p;
M = recovStruct.img_size;
n_f = recovStruct.subframe_l;
upsamplef = recovStruct.upsample_factor;

%number of boundry pixels to guard against artifacts
n_boundry_p = 5;

%% 2- extract x-y channel parameters

%------------------------------------------------------------
%x_channel Fourier transfor of PSFs
FXX = FPSF.FXX;

%gradients
FXXdx = FPSF.FXXdx;
FXXdy = FPSF.FXXdy;

boundxryIndx = [1:n_boundry_p, sqrt(N) - n_boundry_p + 1:sqrt(N)];

%% 3- computations are performed in Fourier domain

%------------------------------------------------------------
c_temp = (down_sample(fast_mul_fft(reshapeMat(gamma))));

c = bsxfun(@plus, reshape(c_temp, M^2, n_f), b);

g = fast_transp_mul_fft_x((1-SMLM_img./c));

%% utility functions

%------------------------------------------------------------


    function out = reshapeMat(A)
        out = reshape(A, sqrt(N), sqrt(N), 3, n_f);
    end

    function down_sampled = down_sample(x)
        down_sampled = x(1:upsamplef:end, 1:upsamplef:end, :);
    end


    function out_N = fast_mul_fft(x)

        out_N = real(ifft2(bsxfun(@times, FXX, fft2(xxgrid(x))) +...
            bsxfun(@times, FXXdx, fft2(xxdxgrid(x))) + bsxfun(@times, FXXdy, fft2(xxdygrid(x)))));
    end


    function out_N1_inN2_t = xxgrid(x)
        out_N1_inN2_t = (reshape(x(:, :, 1, :), sqrt(N), sqrt(N), n_f));
        %         out_N1_inN2=(padarray(out_N1_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N1_inN2_t(boundxryIndx, :) = 0;
        out_N1_inN2_t(:, boundxryIndx) = 0;

    end

    function out_N3_inN2_t = xxdxgrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 2, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2_t=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end
    function out_N3_inN2_t = xxdygrid(x)
        out_N3_inN2_t = (reshape(x(:, :, 3, :), sqrt(N), sqrt(N), n_f));
        %         out_N3_inN2=(padarray(out_N3_inN2_t(6:end-5,6:end-5,:),[5,5]));
        out_N3_inN2_t(boundxryIndx, :) = 0;
        out_N3_inN2_t(:, boundxryIndx) = 0;
    end

    function out_N3 = fast_transp_mul_fft_x(z)
        z_up = up_samp(z);
        fz = fft2(z_up);
        xx_temp = ifft2(bsxfun(@times, conj(FXX), fz));
        xxdx_temp = ifft2(bsxfun(@times, conj(FXXdx), fz));
        xxdy_temp = ifft2(bsxfun(@times, conj(FXXdy), fz));

        indx = [1:n_boundry_p, size(xx_temp, 1) - n_boundry_p + 1:size(xx_temp, 1)];
        xx_temp(indx, :) = 0;
        xx_temp(:, indx) = 0;
        xxdx_temp(:, indx) = 0;
        xxdy_temp(:, indx) = 0;


        out_N3 = real([reshape(xx_temp, N, n_f);
            reshape(xxdx_temp, N, n_f); ...
            reshape(xxdy_temp, N, n_f)]);

        function out_N1_N3 = up_samp(x)
            out_N1_N3 = zeros(sqrt(N), sqrt(N), n_f);
            out_N1_N3(1:upsamplef:end, 1:upsamplef:end, :) = reshape(x, M, M, n_f);
        end


    end

end