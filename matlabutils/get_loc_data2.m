function loc_data = get_loc_data2(gammaf, recovStruct)
%get_loc_data computes molecular parameters (localization data) consist of
%brightness, position and confidence level
%->---
%input
%->---
%gammaf:    array(3N,n_f) -molecular parameter estimates
%grid_p:    array(N,1)    -lateral grid points along x axis
%T:         array(N,N)    -point-wise confidence map
%thr:       scalar [0,1]  -point-wise confidence threshold
%---->-
%output
%---->-
%loc_data:   array with columns [frame _number  brightness, x,y, confidence]

loc_data = [];

%grid points passing through origin along x
grid_p = recovStruct.lateral_grid_p;

%number of grid points along x axis
grid_length = length(grid_p);

%total number of grid points
n_grid_p = recovStruct.n_grid_p;


%number of frames
n_frames = recovStruct.num_frames;

%loop over frames
for i = 1:n_frames

    %current molecular estimates
    z = gammaf(:, i);

    %molecular parameters related to XX basis
    z_xx = [z(1:n_grid_p); z(1 * n_grid_p + 1:2 * n_grid_p); z(2 * n_grid_p + 1:3 * n_grid_p)];

    %extract brightness-scaled second moments estimates (XX, YY and ZZ)
    z_secondMoments = reshape(z(1:n_grid_p), sqrt(n_grid_p), sqrt(n_grid_p), 1);

    %checking if current frame contains molecular parameters
    if any(z_secondMoments(:) > 0)


        % get the index position of molecules across XX, YY and ZZ
        [I, J] = find(sum(z_secondMoments, 3) > 0);
        loc_t = [J, I];


        k = 1;
        while k <= size(loc_t, 1)

            pm = [grid_p(loc_t(k, 1)), grid_p(loc_t(k, 2))]; %continuous location
            im = (loc_t(k, 1) - 1) * grid_length + loc_t(k, 2); %index

            % XX basis
            %-----------------------------------------------------


            imx_xx = pm(1, 1) + (z_xx(im + n_grid_p) / (eps + z_xx(im))) * 10^2; %x position
            imy_xx = pm(1, 2) + (z_xx(im + 2 * n_grid_p) / (eps + z_xx(im))) * 10^2; %y position
            br_m_xx = z_xx(im); %brightness scales

            % combine position estimates
            %-----------------------------------------------------

            br_m = br_m_xx; % molecule brightness is the sum across XX, YY and ZZ

            imx = imx_xx;
            imy = imy_xx;

            % map to sencond moments
            %-----------------------------------------------------


            %update localizaton data
            %----------------------------------------------------
            loc_data = [loc_data; i, imx, imy, br_m];

            k = k + 1;
        end

    end

end