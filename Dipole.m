%%Dipole Model Magnetic Geometry Helper Function for Dipole Model.
function [Bx,By,Bz] = Dipole(x,y,z,season,rot,planet)

Re=1; %Radius of planet (normalized to planetary radii)

if planet == 0
    raxis = 97.8; % Tilt of Uranus' rotation axis
    mag = -60; % Tilt of Neptune's magnetic axiss
    Me= -2.3e-05; % (T) Surface equatorial magnetic field strength of Uranus
elseif planet == 1
    Me= -14*10^-6; % (T) Surface equatorial magnetic field strength of Neptune. negative sign is so the field shoots out the top
    raxis = -28.3; % Tilt of Neptune's rotation axis
    mag = -46.8; % Tilt of Neptune's magnetic axis
end

% 1. rotate for season
if season == 1 % solstice, tilt around y
    xra = x*cosd(raxis) + z*sind(raxis);
    yra = y;
    zra = -x*sind(raxis) + z*cosd(raxis);
elseif season == 0 % equinox, tilt around x
    xra = x;
    yra = y*cosd(raxis) - z*sind(raxis);
    zra = y*sind(raxis) + z*cosd(raxis);
else
    disp('Season is listed incorrectly. Check for error.')
end


% 2. rotate magnetic axis around north pole 
%rot = -180;

xz = xra*cosd(rot) - yra*sind(rot);
yz = xra*sind(rot) + yra*cosd(rot);
zz = zra; 


% 3. retilt by rotation away from z around y

xp = xz*cosd(mag) + zz*sind(mag);
yp = yz;
zp = -xz*sind(mag) + zz*cosd(mag);


% solve for B
Bx1=((-3*Me*xp.*zp.*Re^3))./((sqrt(xp.^2+yp.^2+zp.^2)).^5);
By1=((-3*Me*yp.*zp.*Re^3))./((sqrt(xp.^2+yp.^2+zp.^2)).^5);
Bz1=((Me*(Re^3)*(-2*zp.^2+xp.^2+yp.^2)))./((sqrt(xp.^2+yp.^2+zp.^2)).^5);


% 3.
Bx2 = Bx1*cosd(mag) - Bz1*sind(mag);
By2 = By1;
Bz2 = Bx1*sind(mag) + Bz1*cosd(mag);
% 2.
Bx3 = Bx2*cosd(rot) + By2*sind(rot);
By3 = -Bx2*sind(rot) + By2*cosd(rot); 
Bz3 = Bz2;

if season == 1 % solstice, tilt around y
    Bx4 = Bx3*cosd(raxis) - Bz3*sind(raxis);
    By4 = By3;
    Bz4 = Bx3*sind(raxis) + Bz3*cosd(raxis);
elseif season == 0 % equinox, tilt around x
    Bx4 = Bx3;
    By4 = By3*cosd(raxis) + Bz3*sind(raxis);
    Bz4 = -By3*sind(raxis) + Bz3*cosd(raxis);
else
    disp('Season is listed incorrectly. Check for error.')
end

Bx = Bx4;
By = By4;
Bz = Bz4;

end

