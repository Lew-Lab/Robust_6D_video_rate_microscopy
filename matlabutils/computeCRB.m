function [CRBr,CRBm] = computeCRB(B,Bgradx,Bgrady,Bgradz,s,b,m)

for i = 1:numel(B)
    Btmp = B{i};
    Bscaling = sum(sum(Btmp(:,:,1)));
    Btmp = Btmp/Bscaling;
    img = Btmp(:,:,1)*m(1)+Btmp(:,:,2)*m(2)+Btmp(:,:,3)*m(3)+...
        Btmp(:,:,4)*m(4)+Btmp(:,:,5)*m(5)+Btmp(:,:,6)*m(6);
    s = s/sum(sum(img));
    img = img*s+b;
    img(img<1e-15) = nan;
    Btmp = Bgradx{i}/Bscaling;
    img_dx = (Btmp(:,:,1)*m(1)+Btmp(:,:,2)*m(2)+Btmp(:,:,3)*m(3)+...
        Btmp(:,:,4)*m(4)+Btmp(:,:,5)*m(5)+Btmp(:,:,6)*m(6))*s;
    Btmp = Bgrady{i}/Bscaling;
    img_dy = (Btmp(:,:,1)*m(1)+Btmp(:,:,2)*m(2)+Btmp(:,:,3)*m(3)+...
        Btmp(:,:,4)*m(4)+Btmp(:,:,5)*m(5)+Btmp(:,:,6)*m(6))*s;
    Btmp = Bgradz{i}/Bscaling;
    img_dz = (Btmp(:,:,1)*m(1)+Btmp(:,:,2)*m(2)+Btmp(:,:,3)*m(3)+...
        Btmp(:,:,4)*m(4)+Btmp(:,:,5)*m(5)+Btmp(:,:,6)*m(6))*s;
    FIr(1,1) = nansum(nansum(img_dx.^2./img));
    FIr(2,2) = nansum(nansum(img_dy.^2./img));
    FIr(3,3) = nansum(nansum(img_dz.^2./img));
    FIr(1,2) = nansum(nansum(img_dx.*img_dy./img));
    FIr(1,3) = nansum(nansum(img_dx.*img_dz./img));
    FIr(2,3) = nansum(nansum(img_dz.*img_dy./img));
    FIr(2,1) = FIr(1,2);
    FIr(3,1) = FIr(1,3);
    FIr(3,2) = FIr(2,3);
    CRBr(:,:,i) = inv(FIr);
    Btmp = B{i};
    Btmp = Btmp/Bscaling;
    for j = 1:6
        for k = j:6
            FIm(j,k) = s^2*nansum(nansum(Btmp(:,:,j).*Btmp(:,:,k)./img));
            FIm(k,j) = FIm(j,k);
        end
    end
    CRBm(:,:,i) = inv(FIm);
end

end
