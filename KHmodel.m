% This analytical model evaluates the Kelvin-Helmholtz Instability
% inequality on a magnetopause surface to see if the KHI is possible at
% that current season, rotation, Interplanetary Magnetic Field (IMF)
% strength, and IMF direction.s

clc
close all
clear all

fid = fopen( 'KHModelData.csv', 'w' );
fprintf( fid, '%s,%s,%s,%s,%s,%s,%s\n', 'Planet', 'season', 'k-value', 'IMF_dir', 'IMF_mag', 'Rotation_Angle', 'Percent_Area');

%% Constants
pmass = 1.67*10^-27; % (kg) proton mass
mu_0 = 4*pi*10^-7; % (H/m) vacuum permeability
emass = 9.109 * 10^-31; % (kg) electron mass
q = -1.602*10^-19; % (C) electron charge
flowspeed = 400 *10^3; % (m/s) speed of solar wind
kb = 1.38*10^-23; % (JK-1) Boltzmann's constant
gamma = 5/3; % (dimensionless) used in sonic mach number calculation

%% Define the parameters to change each run.

seasons = [0 1]; % Define Equinox = 0 and Solstice = 1
kval = [0.75 2 200]; % Magnitude of the k-vector
IMFdir = [0 1 2]; % Define IMF_z = 0, IMF_y = 1, and IMF_zy = 2
IMFmags = [0.01e-9 0.1e-9 0.3e-9 0.5e-9 1e-9]; % Set IMF in Tesla
rotdir = [0 45 90 135 180 225 270 315 360]; % The diurnal rotation of the planet
planets=[0 1]; % Define which planet to test: Uranus = 0 and Neptune = 1

for hp = 1:length(planets)
    for ip = 1:length(seasons)
        for jp = 1:length(kval)
            datvar = zeros([length(IMFdir),length(IMFmags),length(rotdir)]);
            for kp = 1:length(IMFdir)
                for lp = 1:length(IMFmags)
                    for mp = 1:length(rotdir)

                        close all % Close any figures open from the previous iteration

                        % IMF direction being used for this iteration
                        if IMFdir(kp) == 0
                            IMF_z = IMFmags(lp);
                            IMF_y = 10e-200;
                            IMF_x = 10e-200;
                        elseif IMFdir(kp) == 1
                            IMF_y = IMFmags(lp);
                            IMF_z = 10e-200;
                            IMF_x = 10e-200;
                        elseif IMFdir(kp) == 2
                            IMF_z = IMFmags(lp)./sqrt(2);
                            IMF_y = IMFmags(lp)./sqrt(2);
                            IMF_x = 10e-200;
                        end

                        if planets(hp) == 0

                            planet = 0;
                            planetstr = 'Uranus';

                            % Uranus specific constants
                            rotperiod = 17.23*3600;
                            R = 25362; % (km) radius of Uranus
                            MP = 20; % magnetopause standoff distance in planetary radii
                            BS = 23.7; % bowshock standoff distance in planetary radii
                            limitplotting = 68; % V2 exited Uranus magnetosphere at 68 RU
                            surface_den = 0.002 * 10^6*pmass; % (kg/m3) from Fig 3 (https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/91JA01857)
                            density_sw = 0.03 * 10^6; % (1/m^3)
                            density_sw = density_sw * pmass; % (kg/m^3)
                            raxis = 97.8; % Tilt of Uranus' rotation axis
                            MS = 27; % Sonic Mach number

                        elseif planets(hp)==1
                            planet =1;
                            planetstr='Neptune';

                            % Neptune specific constants
                            rotperiod = 16.1*3600;
                            R = 24622; % (km) radius of Neptune
                            MP = 26; % magnetopause standoff distance in planetary radii
                            BS = 34.9; % bowshock standoff distance in planetary radii
                            limitplotting = 60; % V2 exited neptune magnetosphere at 60 RN
                            surface_den = 0.01 * (0.5*10^6*pmass + 0.5*10^6*14.0067*pmass); % (kg/m3) from Fig 3 (https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/91JA01857) 50% h+ and 50% N+
                            density_sw = 0.025 * 10^6; % (1/m^3)
                            density_sw = density_sw * pmass; % (kg/m^3)
                            raxis = -28.3; % Tilt of Neptune's rotation axis
                            MS = 28; % Sonic Mach number
                        end

                        % Wave number k being used.
                        kmagnitude = kval(jp);

                        % Season for this iteration
                        season = seasons(ip); 

                        %%  Set up volume
                        % all grid spaces are in units of one planetary radii
                        n = 201; % resolution
                        xi = linspace(-MP-1,limitplotting,n); % max and min are determined by looking at max and min computed by using yi and zi in paraboloid equation below
                        yi = linspace(-70,70,n);
                        zi = linspace(-70,70,n);
                        [xg,yg,zg] = meshgrid(xi,yi,zi);

                        %% Create paraboloid as magnetopause surface
                        [yii,zii] = meshgrid(yi,zi);
                        paraboloid = 1/(2*MP)*(yii.^2 + zii.^2) - MP;

                        % Create paraboloids on either side of "magnetopause" to map the 3D grid
                        % values for solar wind and magnetosphere to. This approach is necessary so
                        % the curved grid spaces line up with each other in each paraboloid.
                        mp1 = MP - 0.3;
                        mp2 = MP + 0.3;
                        p1 = 1/(2*mp1)*(yii.^2 + zii.^2) - mp1;
                        p2 = 1/(2*mp2)*(yii.^2 + zii.^2) - mp2;

                        yp = yii;
                        zp = zii;

                        % Set grid at X =limitplotting to be white and therefore invisible.
                        color = zeros(n,n);
                        color1 = zeros(n,n);
                        color2 = zeros(n,n);
                        for j = 1:n
                            for k = 1:n
                                if paraboloid(j,k) > limitplotting
                                    paraboloid(j,k) = limitplotting;
                                    p1(j,k) = limitplotting;
                                    p2(j,k) = limitplotting;
                                    color(j,k) = 0;
                                    color1(j,k) = 0;
                                    color2(j,k) = 0;
                                end

                                if paraboloid(j,k) < limitplotting
                                    color(j,k) = 1;
                                    color1(j,k) = 1;
                                    color2(j,k) = 1;
                                end

                            end
                        end

                        %% Find normal and tangential directions on paraboloid surface. K vector

                        % The normal to the surface
                        [nx,ny,nz] = surfnorm(paraboloid,yp,zp);

                        % Create vectors pointing in negative x direction to use in dot product for
                        % calculating psi angle.
                        xx = -1 * ones(n,n);
                        xy = zeros(n,n);
                        xz = zeros(n,n);

                        % Slope aka tangent on surface gives k vector direction
                        [sy,sz] = gradient(paraboloid);
                        ky1 = sy/max(max(sy));
                        kz1 = sz/max(max(sz));
                        kx1 = zeros(n,n);
                        for j = 1:n
                            for k = 1:n
                                kx1(j,k) = abs(sqrt( 1 - ky1(j,k)^2 - kz1(j,k)^2)); % take absolute value because will get complex number otherwise
                            end
                        end

                        % Multiply k vector components by wave number magnitude
                        kx = kx1 * kmagnitude ;
                        ky = ky1 * kmagnitude;
                        kz = kz1 * kmagnitude;

                        %% Define grid spaces where solar wind and magnetosphere are located.
                        grid = zeros(n,n,n);
                        for i = 1:n
                            for j = 1:n
                                for k = 1:n
                                    if xg(1,j,1) < paraboloid(i,k)
                                        grid(i,j,k) = 1; % solar wind space
                                    else
                                        grid(i,j,k) = 0; % magnetosphere space
                                    end
                                end
                            end
                        end


                        %% Designate density and magnetic field for magnetosphere and solar wind
                        % Initialize vectors
                        swb_z = zeros(n,n,n);
                        swb_x = zeros(n,n,n);
                        swb_y = zeros(n,n,n);
                        sw_den = zeros(n,n,n);
                        l = zeros(n,n,n);
                        drape_bx = zeros(n,n,n);
                        drape_by = zeros(n,n,n);
                        drape_bz = zeros(n,n,n);
                        drape_tot = zeros(n,n,n);
                        swb_tot = zeros(n,n,n);
                        planet_bx = zeros(n,n,n);
                        planet_by = zeros(n,n,n);
                        planet_bz = zeros(n,n,n);
                        planet_den = zeros(n,n,n);

                        % Set initial solar wind values in 3D space
                        for i = 1:1:n
                            for j = 1:1:n
                                for k = 1:1:n
                                    if grid(i,j,k) == 0 % if magnetosphere, set solar wind values to zero at these spots in grid
                                        swb_z(i,j,k) = 0;
                                        swb_x(i,j,k) = 0;
                                        swb_y(i,j,k) = 0;
                                        sw_den(i,j,k) = 0;
                                    else % if solar wind, set solar wind values in grid
                                        swb_z(i,j,k) = IMF_z;
                                        swb_y(i,j,k) = IMF_y;
                                        swb_x(i,j,k) = IMF_x;
                                        swb_tot(i,j,k) = sqrt(swb_z(i,j,k)^2 + swb_y(i,j,k)^2 + swb_x(i,j,k)^2);
                                        sw_den(i,j,k) = density_sw;
                                    end
                                end
                            end
                        end

                        % Employ approach from Cooling et al (2001) to drape the IMF field lines
                        A = (2*BS - MP)/(2*(BS-MP));
                        for i = 1:1:n
                            for j = 1:1:n
                                for k = 1:1:n
                                    l(i,j,k) = xg(i,j,k);
                                    drape_bx(i,j,k) = -A*(-swb_x(i,j,k)*(1-MP/(2*l(i,j,k))) + swb_y(i,j,k)*-1*yg(i,j,k)/l(i,j,k) + swb_z(i,j,k)*zg(i,j,k)/l(i,j,k) );
                                    drape_by(i,j,k) =  A*(-swb_x(i,j,k)*-1*yg(i,j,k)/(2*l(i,j,k)) + swb_y(i,j,k)*(2-(-1*yg(i,j,k))^2/(l(i,j,k)*MP) ) - swb_z(i,j,k)*-1*yg(i,j,k)*zg(i,j,k)/(l(i,j,k)*MP));
                                    drape_bz(i,j,k) = A*(-swb_x(i,j,k)*zg(i,j,k)/(2*l(i,j,k)) + swb_y(i,j,k)*-1*yg(i,j,k)*zg(i,j,k)/(l(i,j,k)*MP) + swb_z(i,j,k)*(2 - (zg(i,j,k))^2/(l(i,j,k)*MP)));

                                    % Modify for space where x > 0
                                    if xg(i,j,k) > 0
                                        drape_bz(i,j,k) = -1*drape_bz(i,j,k);
                                        drape_bx(i,j,k) = -1*drape_bx(i,j,k);
                                        drape_by(i,j,k) = -1*drape_by(i,j,k);
                                    end

                                    drape_tot(i,j,k) = sqrt(drape_bx(i,j,k)^2 + drape_by(i,j,k)^2 + drape_bz(i,j,k)^2);

                                    % catch for if the draped field total is somehow larger than
                                    % the original IMF
                                    if drape_tot(i,j,k) > swb_tot(i,j,k)
                                        drape_bx(i,j,k) = (drape_bx(i,j,k)/drape_tot(i,j,k)) * swb_tot(i,j,k);
                                        drape_by(i,j,k) = (drape_by(i,j,k)/drape_tot(i,j,k)) * swb_tot(i,j,k);
                                        drape_bz(i,j,k) = (drape_bz(i,j,k)/drape_tot(i,j,k)) * swb_tot(i,j,k);
                                    end


                                end
                            end
                        end


                        % Set initial solar wind values in 3D space
                        % make dipole everywhere by calling dipole function script Joe Cags created
                        % the tilt is included here
                        [planetbx,planetby,planetbz]= Dipole(xg,yg,zg,season,rotdir(mp),planet);

                        % Set magnetosphere values in 3D space
                        for i = 1:1:n
                            for j = 1:1:n
                                for k = 1:1:n

                                    if grid(i,j,k) == 0 % if magnetosphere, set magnetosphere values in grid

                                        % density
                                        planet_den(i,j,k) = surface_den; % density set as what V2 read immediately after crossing MP for electrons

                                        % magnetic field
                                        planet_bx(i,j,k) = planetbx(i,j,k);
                                        planet_by(i,j,k) = planetby(i,j,k);
                                        planet_bz(i,j,k) = planetbz(i,j,k);


                                    else % if outside of magnetosphere, set magnetosphere values to zero at these spots in grid
                                        planet_bx(i,j,k) = 0;
                                        planet_by(i,j,k) = 0;
                                        planet_bz(i,j,k) = 0;

                                    end
                                end
                            end
                        end

                        % Add solar wind values to magnetosphere values to get total field for
                        % plotting
                        btot_x = drape_bx + planet_bx;
                        btot_y = drape_by + planet_by;
                        btot_z =  drape_bz + planet_bz;
                        btot = sqrt( btot_x.^2 + btot_y.^2 + btot_z.^2 );
                        btot = log10(btot); % differences more obvious when plotting

                        %% Set Magnetosphere Plasma Velocity

                        % initialize vectors
                        mag_v = zeros(n,n,n);
                        r = zeros(n,n,n);
                        xgs = zeros(n,n,n);
                        ygs = zeros(n,n,n);
                        zgs = zeros(n,n,n);

                        t = rotperiod; % (s) Planet rotation period
                        w = 2*pi/t; % (rad/s) angular velocity

                        for i = 1:n
                            for j = 1:n
                                for k = 1:n

                                    % tilt based on season
                                    if season == 1 % solstice
                                        xgs(i,j,k) = xg(i,j,k)*cosd(raxis) + zgs(i,j,k)*sind(raxis);
                                        ygs(i,j,k) = yg(i,j,k);
                                        zgs(i,j,k) = -xg(i,j,k)*sind(raxis) + zg(i,j,k)*cosd(raxis);
                                    end

                                    if season == 0 % equinox
                                        xgs(i,j,k) = xg(i,j,k);
                                        ygs(i,j,k) = yg(i,j,k)*cosd(raxis) - zg(i,j,k)*sind(raxis);
                                        zgs(i,j,k) = ygs(i,j,k)*sind(raxis) + zgs(i,j,k)*cosd(raxis);
                                    end

                                    % compute rigid corotational velocity
                                    r(i,j,k) = sqrt(xgs(i,j,k)^2 + ygs(i,j,k)^2 + zgs(i,j,k)^2);
                                    if grid(i,j,k) == 0 % if magnetosphere
                                        mag_v(i,j,k) = w*r(i,j,k)*R*10^3; % (m/s) because each grid point is 1 RN
                                    end

                                    % reduce corotational plasma based on Cassini data from Wilson et al. (2009)
                                    mag_v(i,j,k) = 0.8 * mag_v(i,j,k);

                                end
                            end
                        end

                        %% Solar Wind Plasma Velocity
                        % calculation based on equation from Masters (2014)
                        % initialize vectors
                        dotsum = zeros(n,n);
                        magnorm = zeros(n,n);
                        psi = zeros(n,n);
                        sw_v = zeros(n,n);

                        % calculate psi angle
                        for i = 1:n
                            for j = 1:n
                                dotsum(i,j) = nx(i,j)*xx(i,j) + ny(i,j)*xy(i,j) + nz(i,j)*xz(i,j);
                                magnorm(i,j) = (nx(i,j)^2 + ny(i,j)^2 + nz(i,j)^2)^0.5;
                                psi(i,j) = acosd( dotsum(i,j)/magnorm(i,j) );
                                % Take the difference between psi and 180 to define the x axis in the correct direction
                                psi(i,j) = 180 - psi(i,j); % psi is angle with x in negative direction and the normal vector

                                sw_v(i,j) = flowspeed * sqrt( (MS^2 +3)/MS^2 * ( 1 - ( (cosd(psi(i,j)))^2 + 3^(5/2)/4^4 * (5*MS^2-1)^(3/2)/MS^5 * sind(psi(i,j))^2 )^(2/5) ) );

                            end
                        end

                        % Produce Figure 3
                        % figure
                        % imagesc(sw_v)

                        %% Assign 3D density and magnetic field values to the 2D paraboloid surface

                        % initialize matrices only using values that lie right along the magnetopause
                        % border
                        bsw_den = zeros(n,n);
                        bswb_z = zeros(n,n);
                        bswb_y = zeros(n,n);
                        bswb_x = zeros(n,n);
                        bswb_tot = zeros(n,n);
                        bplanet_den = zeros(n,n);
                        bplanet_bx = zeros(n,n);
                        bplanet_by = zeros(n,n);
                        bplanet_bz = zeros(n,n);
                        test_p1 = zeros(n,n,n);
                        test_p2 = zeros(n,n,n);
                        bswb_zd = zeros(n,n);
                        bswb_yd = zeros(n,n);
                        bswb_xd = zeros(n,n);
                        bswb_totd = zeros(n,n);

                        for i = 1:n
                            for j = 1:n
                                for k = 1:n

                                    % difference between position in x direction and paraboloids
                                    % bordering magnetopause
                                    test_p1(i,j,k) = abs(xg(1,j,1) - p1(i,k));
                                    test_p2(i,j,k) = abs(xg(1,j,1) - p2(i,k));


                                    difference = 1.50;

                                    if grid(i,j,k) == 0 % magnetosphere space
                                        bplanet_den(i,k) = surface_den; % constant in space inside the magnetosphere
                                        if test_p1(i,j,k) < difference
                                            bplanet_bx(i,k) = planet_bx(i,j,k);
                                            bplanet_by(i,k) = planet_by(i,j,k);
                                            bplanet_bz(i,k) = planet_bz(i,j,k);
                                        end

                                    elseif grid(i,j,k) == 1 % solar wind space
                                        bsw_den(i,k) = density_sw;
                                        if test_p2(i,j,k) < difference
                                            bswb_zd(i,k) = drape_bz(i,j,k);
                                            bswb_yd(i,k) = drape_by(i,j,k);
                                            bswb_xd(i,k) = drape_bx(i,j,k);
                                            bswb_totd(i,k) = sqrt(bswb_xd(i,k)^2 + bswb_yd(i,k)^2 + bswb_zd(i,k)^2);
                                        end

                                        % adjust solar wind density and IMF for the PDL
                                        bsw_den(i,k) = sqrt(0.15) * bsw_den(i,k);
                                        bswb_z(i,k) = (0.15)^(-0.5)* bswb_zd(i,k);
                                        bswb_y(i,k) = (0.15)^(-0.5)* bswb_yd(i,k);
                                        bswb_x(i,k) = (0.15)^(-0.5)* bswb_xd(i,k);
                                        bswb_tot(i,k) = sqrt(bswb_x(i,k)^2 + bswb_y(i,k)^2 + bswb_z(i,k)^2);

                                    end

                                end

                            end
                        end

                        %% Calculate the right hand side of the condition
                        % initialize
                        kdotb1 = zeros(n,n);
                        kdotb2 = zeros(n,n);
                        addkdotb = zeros(n,n);
                        rightside = zeros(n,n);

                        for j = 1:n
                            for k = 1:n
                                kdotb1(j,k) = kx(j,k)*bswb_x(j,k) + ky(j,k)*bswb_y(j,k) + kz(j,k)*bswb_z(j,k);
                                kdotb2(j,k) = kx(j,k)*bplanet_bx(j,k) + ky(j,k)*bplanet_by(j,k) + kz(j,k)*bplanet_bz(j,k);
                                addkdotb(j,k) = kdotb1(j,k)^2 + kdotb2(j,k)^2;

                                rightside(j,k) = 1/mu_0 * (1/bsw_den(j,k) + 1/bplanet_den(j,k) ) * addkdotb(j,k);
                            end
                        end

                        rightside(isinf(rightside)) = 10e20; % set any values that are infinity to this high number
                        
                        % Produce plots in Figure 8
                        % figure
                        % imagesc(kdotb2)

                        %% Assign 3D velocity values to 2D paraboloid surface
                        % initialize matrices
                        mag_v0 = zeros(n,n);
                        vdif = zeros(n,n);
                        mv_theta = zeros(n,n);
                        mv_psi = zeros(n,n);
                        mag_vx = zeros(n,n);
                        mag_vy = zeros(n,n);
                        mag_vz = zeros(n,n);
                        modified_p = zeros(n,n);
                        swtot = zeros(n,n);
                        sw_vx = zeros(n,n);
                        sw_vy = zeros(n,n);
                        sw_vz = zeros(n,n);
                        sw_vtot = zeros(n,n);

                        % overall magnetosphere velocity to 2D
                        for i = 1:n
                            for j = 1:n
                                for k = 1:n
                                    if grid(i,j,k) == 0 % magnetosphere space
                                        if test_p1(i,j,k) < difference % magnetosphere parabola, diffference variable is defined above
                                            mag_v0(i,k) = mag_v(i,j,k);
                                        end
                                    end
                                end
                            end
                        end

                        % now define magnetosphere velocity components for 2D surface points
                        for j = 1:n
                            for k = 1:n

                                mv_theta(j,k) = acosd(zp(j,k)/sqrt(paraboloid(j,k)^2 + yp(j,k)^2 + zp(j,k)^2));
                                if abs(paraboloid(j,k)) < 0.5 && yp(j,k) > 0
                                    mv_psi(j,k) = 90;
                                elseif abs(paraboloid(j,k)) < 0.5 && yp(j,k) < 0
                                    mv_psi(j,k) = -90;
                                elseif paraboloid(j,k) > 0
                                    mv_psi(j,k) = atand(yp(j,k)/paraboloid(j,k));
                                elseif paraboloid(j,k) < 0
                                    mv_psi(j,k) = atand(yp(j,k)/paraboloid(j,k)) + 180;
                                else
                                end


                                % use angles to calculate how much of magnetopause velocity goes to
                                % each component
                                mag_vx(j,k) = mag_v0(j,k) * cosd(mv_psi(j,k)) * sind(mv_theta(j,k));
                                mag_vy(j,k) = mag_v0(j,k) * sind(mv_psi(j,k)) * sind(mv_theta(j,k));
                                mag_vz(j,k) = mag_v0(j,k) * cosd(mv_theta(j,k));


                            end
                        end

                        % define solar wind velocity components for 2D surface points
                        for j = 1:n
                            for k = 1:n
                                % use same angles to calculate how much of solar wind velocity goes
                                % to each directional component
                                sw_vy(j,k) = sw_v(j,k) * sind(mv_psi(j,k)) * sind(mv_theta(j,k));
                                sw_vz(j,k) = sw_v(j,k) * cosd(mv_theta(j,k));
                                sw_vx(j,k) = sw_v(j,k)* cosd(mv_psi(j,k)) * sind(mv_theta(j,k));

                                sw_vtot(j,k) = sqrt(sw_vx(j,k)^2 + sw_vy(j,k)^2 + sw_vz(j,k)^2);

                            end
                        end

                        %% Calculate the left side of the equation
                        % initialize
                        vdiff = zeros(n,n);
                        leftside = zeros(n,n);

                        for i = 1:n
                            for k = 1:n
                                vdiff(i,k) = sqrt((sw_vx(i,k) - mag_vx(i,k))^2 + (sw_vy(i,k) - mag_vy(i,k))^2 + (sw_vz(i,k) - mag_vz(i,k))^2);
                                leftside(i,k) = ( kx(i,k)*(sw_vx(i,k) - mag_vx(i,k)) + ky(i,k)*(sw_vy(i,k) - mag_vy(i,k)) + kz(i,k)*(sw_vz(i,k) - mag_vz(i,k)) )^2;

                            end
                        end


                        %% Final result and plot
                        % initialize
                        total = 0;
                        true = 0;
                        plotresult = zeros(n,n);

                        for j = 1:n
                            for k = 1:n
                                if paraboloid(j,k) < limitplotting
                                    total = total +1;

                                    if leftside(j,k) >= rightside(j,k)
                                        true = true +1;
                                        plotresult(j,k) = 1;
                                    end

                                    if rightside(j,k) > leftside(j,k)
                                        plotresult(j,k) = 0;
                                    end

                                end

                                if paraboloid(j,k) >= limitplotting
                                    plotresult(j,k) = 2;

                                end

                            end
                        end


                        % display surface area percentage in command window
                        percent = (true./total).*100;
                        datvar(kp,lp,mp) = percent;

                        cmap1 = [0 0 0
                            1 0 1
                            1 1 1];
                        cmap2 =[0 0 0
                            1 0 1
                            1 1 1];

                        figure()
                        surf(paraboloid,yp,zp,plotresult,'EdgeColor','none')
                        colormap(cmap1)
                        title('Regions where KH condition is satisfied')
                        axis equal
                        view(-90,0) % look down YZ plane

                        if season == 1
                            seasonstr = '_Solstice';
                            seasonstr2 = 'Solstice';
                        elseif season == 0
                            seasonstr = '_Equinox';
                            seasonstr2 = 'Equinox';
                        end

                        if IMF_z > IMF_y
                            IMFstr = '_IMF_z_';
                            IMFstr2 = 'z';
                        elseif IMF_y > IMF_z
                            IMFstr = '_IMF_y_';
                            IMFstr2 = 'y';
                        elseif IMF_y == IMF_z
                            IMFstr = '_IMF_zy_';
                            IMFstr2 = 'zy';
                        end
                        IMFmag = sqrt(IMF_x.^2+IMF_y.^2+IMF_z.^2).*10^9;
                        set(gcf, 'visible', 'off');
                        fprintf( fid, '%s,%s,%d,%s,%f,%d,%f\n', planetstr, seasonstr2, kmagnitude, IMFstr2, IMFmag, rotdir(mp), percent);
                        print(gcf, [planetstr '_k=' num2str(kmagnitude) seasonstr IMFstr 'IMFmag=' num2str(IMFmag) 'RotAngle=' num2str(rotdir(mp)) 'Area%=' num2str(percent) '.png'],'-dpng','-r600')

                    end

                end

            end

%% Plot the total percent of the surface area where KHIs are possible

            close all

%             red = [0.6350 0.0780 0.1840];
%             blue = [0 0.4470 0.7410];
%             black = [0 0 0];
%             yellow = [0.9290 0.6940 0.1250];
%             green = [0.4660 0.6740 0.1880];

            colors = {[0.6350 0.0780 0.1840] [0 0.4470 0.7410] [0 0 0] [0.9290 0.6940 0.1250] [0.4660 0.6740 0.1880]};

            for lp = 1:length(IMFmags)
                plot(rotdir,squeeze(datvar(1,lp,:)),'LineWidth',2.0,'color',colors{lp},'DisplayName',['|B| = ' num2str(IMFmags(lp).*10.^9) ' nT'])
                hold on
                plot(rotdir,squeeze(datvar(2,lp,:)),':','LineWidth',2.0,'color',colors{lp},'HandleVisibility','off')
                plot(rotdir,squeeze(datvar(3,lp,:)),'--','LineWidth',2.0,'color',colors{lp},'HandleVisibility','off')
            end


            hold off
            set(gcf, 'visible', 'off'); % set the 'visible' property of the figure to 'off'


            leg = legend;
            leg.Color = 'k';
            ylim([0 100])
            xlim([0 360])
            x0=20; %Image offset from corner of screen
            y0=20; %Image offset from corner of screen
            width=750; %Image Resolution (may need to be reduced for smaller screens)
            height=900; %Image Resolution (may need to be reduced for smaller screens)
            x = [x0 y0 width-20 height-20];
            ax = gca;
            ax.TitleFontSizeMultiplier = 1.1;
            set(gcf,'position',[x0 y0 width height]);
            title([planetstr 'at ' seasonstr2 ': KHI Surface Area possible for different IMF. k = ' num2str(kmagnitude)],'color','k');
            xlabel('Rotation','color','k');
            xticks(rotdir);
            ylabel('Surface Area Percent','color','k');
            print(gcf, ['Magnetic_Field_Rotation_Plot_' planetstr '_k=' num2str(kmagnitude) seasonstr '.png'],'-dpng','-r600');
            close all
        end
    end
end

fclose( fid );
   