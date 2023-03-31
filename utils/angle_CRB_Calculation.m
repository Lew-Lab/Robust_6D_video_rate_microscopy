function sigmaDelta = angle_CRB_Calculation(basis_matrix,theta,phi,gamma,signal,background)
% 06/16/2022 Yiyang modified based on Oumeng's computeCRB.m

dtheta = 10^-6;
dphi = 10^-6;

M_orig = angle2M_2(theta,phi,gamma);
M_dtheta = angle2M_2(theta+dtheta,phi,gamma);
M_dphi = angle2M_2(theta,phi+dphi,gamma);

I = signal*basis_matrix*M_orig+background;
I_dtheta = signal*basis_matrix*M_dtheta+background;
I_dphi = signal*basis_matrix*M_dphi+background;

I_grad_theta = (I_dtheta - I)/dtheta;
I_grad_phi = (I_dphi - I)/dphi;

angleFI = zeros(2);

angleFI(1,1) = sum(I_grad_theta.^2./I);
angleFI(2,2) = sum(I_grad_phi.^2./I);
angleFI(1,2) = sum(I_grad_theta.*I_grad_phi./I);
angleFI(2,1) = sum(I_grad_phi.*I_grad_theta./I);

angleCRB = inv(angleFI);
sigmaDelta = 2*asin(sqrt(sin(theta)/(4*pi)*sqrt(det(angleCRB))));

end
