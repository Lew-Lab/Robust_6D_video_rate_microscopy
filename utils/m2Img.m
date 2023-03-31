function I = m2Img(basisMatrix,M,signal,background)
% basisMatrix is a N*6 matrix, representing the 6 basis images
% M is a 6*1 vector, representing the second moments

I = signal.*basisMatrix*M+background;
% yield a N*1 vector, representing the image
end