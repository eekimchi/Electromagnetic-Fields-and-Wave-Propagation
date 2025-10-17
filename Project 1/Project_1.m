clc;
clear;
close all;

% Every plot is technically correct, but everything is off by a phase
%  of some sort by roughly 10 degrees

% conductivity
sigma = 0;  % conductivity is 1 for all materials used

% permeability
mu0 = (4*pi)*(10^-7);
mu = mu0; % permeability is 1 for all materials used

% permittivity
epsilon0 = 8.85418e-12;
epsilon_air = epsilon0;
epsilon_teflon = epsilon0*2.1;
epsilon_custom1 = epsilon0*2;
epsilon_custom2 = epsilon0*3;
epsilon_custom3 = epsilon0*4;
epsilon_custom4 = epsilon0*3.5;
epsilon_custom5 = epsilon0*5;
epsilon_custom6 = epsilon0*6.5;
epsilon_custom7 = epsilon0*8;

% impedance for each region
eta_air = sqrt(mu/epsilon_air);
eta_teflon = eta_air*sqrt(epsilon_air/epsilon_teflon);
eta_custom1 = eta_air*sqrt(epsilon_air/epsilon_custom1);
eta_custom2 = eta_air*sqrt(epsilon_air/epsilon_custom2);
eta_custom3 = eta_air*sqrt(epsilon_air/epsilon_custom3);
eta_custom4 = eta_air*sqrt(epsilon_air/epsilon_custom4);
eta_custom5 = eta_air*sqrt(epsilon_air/epsilon_custom5);
eta_custom6 = eta_air*sqrt(epsilon_air/epsilon_custom6);
eta_custom7 = eta_air*sqrt(epsilon_air/epsilon_custom7);

% length of regions
%   for 1 region, the length is technically infinite
l_base = 0.2;

% frequency
hz = 1e9:0.01e9:5e9;
centerfrequency = 2.5e9;

% convert into radians
omega = hz.*2.*pi;
centeromega = centerfrequency.*2.*pi;

% oblique incidence theta sweep
theta = 0:0.5:90;

% tenth of a wavelength
wavelength10 = ((3e8)/(2.5e9))/10;

% ========== BASE CASE ========== GOOD

% phase constants
%   beta for free space (air)
beta_air = omega.*sqrt(mu.*epsilon_air);

% reflection and transmission
gammaboundary_base = (eta_air - eta_air)/(eta_air + eta_air);
tauboundary_base = (2.*eta_air)/(eta_air + eta_air);

gamma_base = gammaboundary_base.*exp(-1j.*2.*beta_air.*l_base);
tau_base = tauboundary_base.*exp(-1j.*(beta_air.*l_base + beta_air.*l_base));

% Import HFSS results

hfss1 = readmatrix("S Parameter Plot - mag S11 and S21 for Base Case.csv");

figure(1)
plot(hz, abs(gamma_base), hz, abs(tau_base), hz, hfss1(:,2), hz, hfss1(:,3))
xlabel('Frequency in Hz')
ylabel('mag')
ylim([0, 1])
legend('S11','S21','S11 HFSS','S21 HFSS')
title('S Parameter Plot - mag S11 and S21 for Base Case')

hfss2 = readmatrix("S Parameter Plot - deg S11 and S22 Angle Comparison for Base Case.csv");

figure(2)
plot(hz, rad2deg(angle(gamma_base)), hz, rad2deg(angle(gamma_base)), hz, hfss2(:,2), hz, hfss2(:,3))
xlabel('Frequency in Hz')
ylabel('Angle in degrees')
legend('S11','S22','S11 HFSS','S22 HFSS')
title('S Parameter Plot - deg S11 and S22 Angle Comparison for Base Case')


% ========== TWO REGIONS ========== 

% lengths
l_tworegion = 0.1;

% phase constants
beta_custom1 = (omega.*sqrt(mu.*epsilon_custom1));
beta_custom3 = (omega.*sqrt(mu.*epsilon_custom3));
beta_custom7 = (omega.*sqrt(mu.*epsilon_custom7));

% === VARYING PERMITTIVITY === 

% reflection and transmission
gammaboundary_tworegion1 = (eta_custom1 - eta_air)/(eta_custom1 + eta_air);
tauboundary_tworegion1 = ((2*eta_custom1)/(eta_air + eta_custom1))*sqrt(eta_air/eta_custom1);
gammaboundary_tworegion2 = (eta_custom3 - eta_air)/(eta_custom3 + eta_air);
tauboundary_tworegion2 = ((2*eta_custom3)/(eta_air + eta_custom3))*sqrt(eta_air/eta_custom3);
gammaboundary_tworegion3 = (eta_custom7 - eta_air)/(eta_custom7 + eta_air);
tauboundary_tworegion3 = ((2*eta_custom7)/(eta_air + eta_custom7))*sqrt(eta_air/eta_custom7);

gamma_tworegion1 = gammaboundary_tworegion1.*exp(-1j.*2.*beta_air.*l_tworegion);
tau_tworegion1 = tauboundary_tworegion1.*exp(-1j.*(beta_custom1.*l_tworegion + beta_air.*l_tworegion));
gamma_tworegion2 = gammaboundary_tworegion2.*exp(-1j.*2.*beta_air.*l_tworegion);
tau_tworegion2 = tauboundary_tworegion2.*exp(-1j.*(beta_custom3.*l_tworegion + beta_air.*l_tworegion));
gamma_tworegion3 = gammaboundary_tworegion3.*exp(-1j.*2.*beta_air.*l_tworegion);
tau_tworegion3 = tauboundary_tworegion3.*exp(-1j.*(beta_custom7.*l_tworegion + beta_air.*l_tworegion));

% Import HFSS results

hfss3 = readmatrix("S Parameter Plot - mag S11 and S12 for Varying Permittivity for 2 Regions.csv");

figure(3)
plot(hz,abs(gamma_tworegion1),hz,abs(tau_tworegion1),hz,abs(gamma_tworegion2),hz,abs(tau_tworegion2),hz,abs(gamma_tworegion3),hz,abs(tau_tworegion3),hz,hfss3(:,2),hz,hfss3(:,11),hz,hfss3(:,4),hz,hfss3(:,13),hz,hfss3(:,8),hz,hfss3(:,17))
xlabel('Frequency in Hz')
ylabel('mag')
legend('S11 for e=2','S21 for e=2','S11 for e=4','S21 for e=4','S11 for e=8','S21 for e=8', 'S11 for e=2 HFSS','S21 for e=2 HFSS','S11 for e=4 HFSS','S21 for e=4 HFSS','S11 for e=8 HFSS','S21 for e=8 HFSS')
title('S Parameter Plot - mag S11 and S21 for Varying Permittivity for 2 Regions')

% === OBLIQUE INCIDENCE === 

beta_air = centeromega.*sqrt(mu.*epsilon_air);
beta_teflon = centeromega.*sqrt(mu.*epsilon_teflon);

thetat = asind((beta_air.*sind(theta))./beta_teflon);

gammaboundary_tworegion1 = (eta_teflon.*cosd(theta) - eta_air.*cosd(thetat))./(eta_teflon.*cosd(theta) + eta_air.*cosd(thetat));

gamma = gammaboundary_tworegion1.*exp(1j.*2.*beta_air.*l_tworegion.*cosd(theta));
tau = sqrt(1 - abs(gamma).^2);

% Import HFSS results

hfss4 = readmatrix("S Parameter Plot - Perpendicular Polarization Oblique Incidence.csv");

figure(4)
plot(theta, abs(gamma), theta, abs(tau), theta, hfss4(:,2), theta, hfss4(:,3))
xlabel('Incident Angle in Degrees')
ylabel('Magnitude of Coefficients')
legend('S11','S21','S11 HFSS','S21 HFSS')
title('S Parameter Plot - Perpendicular Polarization Oblique Incidence')

gammaboundary_tworegion1 = (-eta_air.*cosd(theta) + eta_teflon.*cosd(thetat))./(eta_air.*cosd(theta) + eta_teflon.*cosd(thetat));

gamma = gammaboundary_tworegion1.*exp(1j.*2.*beta_air.*l_tworegion.*cosd(theta));
tau = sqrt(1 - abs(gamma).^2);

hfss5 = readmatrix("S Parameter Plot - Parallel Polarization Oblique Incidence.csv");

figure(5)
plot(theta, abs(gamma), theta, abs(tau), theta, hfss5(:,2), theta, hfss5(:,3))
xlabel('Incident Angle in Degrees')
ylabel('Magnitude of Coefficients')
legend('S11','S21','S11 HFSS','S21 HFSS')
title('S Parameter Plot - Parallel Polarization Oblique Incidence')

% === VARYING DISTANCE === 

% New varying distances
l_region11 = 4*wavelength10;
l_region12 = 0.200 - l_region11;
l_region21 = 0.1;
l_region22 = 0.200 - l_region21;
l_region31 = 0.131912;
l_region32 = 0.200 - l_region31;

% Beta values
beta_air = omega.*sqrt(mu.*epsilon_air);
beta_teflon = omega.*sqrt(mu.*epsilon_teflon);

% Tau at the Boundary
tauboundary_tworegion1 = ((2*eta_teflon)/(eta_air + eta_teflon)).*sqrt(eta_air/eta_teflon);

tau_tworegion1 = tauboundary_tworegion1.*exp(-1j.*(beta_teflon.*l_region11 + beta_air.*l_region12));
tau_tworegion2 = tauboundary_tworegion1.*exp(-1j.*(beta_teflon.*l_region21 + beta_air.*l_region22));
tau_tworegion3 = tauboundary_tworegion1.*exp(-1j.*(beta_teflon.*l_region31 + beta_air.*l_region32));

% Import HFSS results

hfss6 = readmatrix("S Parameter Plot - S12 Phase for Symmetric vs Asymmetric Regions.csv");

figure(6)
plot(hz, rad2deg(angle(tau_tworegion1)), hz, rad2deg(angle(tau_tworegion2)), hz, rad2deg(angle(tau_tworegion3)), hz, hfss6(:,6), hz, hfss6(:,11), hz, hfss6(:,13))
xlabel('Frequency in Hz')
ylabel('Phase in Degrees')
legend('S21 for l=47.968mm','S21 for l=100mm','S21 for l=131.912mm','S21 for l=47.968mm HFSS','S21 for l=107.928mm HFSS','S21 for l=131.912mm HFSS')
title('S Parameter Plot - S21 Phase for Symmetric vs Asymmetric Regions')

% ========== THREE REGIONS ========== 

% === CHANGING THICKNESS ===

wavelength30 = wavelength10*3;
wavelength50 = wavelength10*5;

beta_custom1 = omega.*sqrt(mu.*epsilon_custom1);

gamma12 = (eta_custom1 - eta_air)/(eta_custom1 + eta_air);
gamma23 = (eta_air - eta_custom1)/(eta_air + eta_custom1);

gamma_tenth = (gamma12 + gamma23.*exp(-1j.*2.*beta_custom1.*wavelength10))./(1 + gamma12.*gamma23.*exp(-1j.*2.*beta_custom1.*wavelength10));
gamma_threetenth = (gamma12 + gamma23.*exp(-1j.*2.*beta_custom1.*wavelength30))./(1 + gamma12.*gamma23.*exp(-1j.*2.*beta_custom1.*wavelength30));
gamma_fivetenth = (gamma12 + gamma23.*exp(-1j.*2.*beta_custom1.*wavelength50))./(1 + gamma12.*gamma23.*exp(-1j.*2.*beta_custom1.*wavelength50));

% Import HFSS results

hfss7 = readmatrix("S Parameter Plot - dB S11 and S12 Varying Thickness for Single Slab.csv");

figure(7)
plot(hz, mag2db(abs(gamma_tenth)), hz, mag2db(abs(gamma_threetenth)), hz, mag2db(abs(gamma_fivetenth)), hz, hfss7(:,2), hz, hfss7(:,4), hz, hfss7(:,6))
xlabel('Frequency in Hz')
ylabel('dB')
legend('S11 for d=12mm','S11 for d=36mm','S11 for d=60mm','S11 for d=12mm HFSS','S11 for d=36mm HFSS','S11 for d=60mm HFSS')
title('S Parameter Plot - dB S11 and S21 Varying Thickness for Single Slab')

% === CHANGING PERMITTIVITY === 

hz = 1e9:0.01e9:5e9;
omega = hz.*2.*pi;

beta_custom1 = omega.*sqrt(mu.*epsilon_custom1);

gamma12 = (eta_custom1 - eta_air)/(eta_custom1 + eta_air);
gamma23 = (eta_air - eta_custom1)/(eta_air + eta_custom1);

gamma_custom1 = (gamma12 + gamma23.*exp(-1j.*2.*beta_custom1.*wavelength10))./(1 + gamma12.*gamma23.*exp(-1j.*2.*beta_custom1.*wavelength10));
tau_custom1 = sqrt(1 - abs(gamma_custom1).^2);

beta_custom3 = omega.*sqrt(mu.*epsilon_custom3);

gamma12 = (eta_custom3 - eta_air)/(eta_custom3 + eta_air);
gamma23 = (eta_air - eta_custom3)/(eta_air + eta_custom3);

gamma_custom3 = (gamma12 + gamma23.*exp(-1j.*2.*beta_custom3.*wavelength10))./(1 + gamma12.*gamma23.*exp(-1j.*2.*beta_custom3.*wavelength10));
tau_custom3 = sqrt(1 - abs(gamma_custom3).^2);

beta_custom5 = omega.*sqrt(mu.*epsilon_custom5);

gamma12 = (eta_custom5 - eta_air)/(eta_custom5 + eta_air);
gamma23 = (eta_air - eta_custom5)/(eta_air + eta_custom5);

gamma_custom5 = (gamma12 + gamma23.*exp(-1j.*2.*beta_custom5.*wavelength10))./(1 + gamma12.*gamma23.*exp(-1j.*2.*beta_custom5.*wavelength10));
tau_custom5 = sqrt(1 - abs(gamma_custom5).^2);

% Import HFSS results

hfss8 = readmatrix("S Parameter Plot - mag S11 and S21 Varying Permittivity for Single Slab.csv");

figure(8)
plot(hz, (abs(gamma_custom1)), hz, (abs(tau_custom1)), hz, (abs(gamma_custom3)), hz, (abs(tau_custom3)), hz, (abs(gamma_custom5)), hz, (abs(tau_custom5)), hz, hfss8(:,2), hz, hfss8(:,7), hz, hfss8(:,4), hz, hfss8(:,9), hz, hfss8(:,5), hz, hfss8(:,10))
xlabel('Frequency in Hz')
ylabel('mag')
legend('S11 for e=2','S21 for e=2','S11 for e=4','S21 for e=4','S11 for e=5','S21 for e=5','S11 for e=2 HFSS','S21 for e=2 HFSS','S11 for e=4 HFSS','S21 for e=4 HFSS','S11 for e=5 HFSS','S21 for e=5 HFSS')
title('S Parameter Plot - mag S11 and S21 Varying Permittivity for Single Slab')

% Changing position doesn't mean anything since in our boundary
%  reflection values, the only beta value we care about is our dielectric

% ========== FIVE REGIONS ==========

% === Air - 3 Dielectrics - Air ===

% Betas
beta_custom1 = omega.*sqrt(mu.*epsilon_custom1);
beta_custom2 = omega.*sqrt(mu.*epsilon_custom2);
beta_custom3 = omega.*sqrt(mu.*epsilon_custom3);
beta_custom7 = omega.*sqrt(mu.*epsilon_custom7);

% === 2-3-4 ===
gamma0 = (eta_custom1 - eta_air)/(eta_custom1 + eta_air);
gamma1 = (eta_custom2 - eta_custom1)/(eta_custom2 + eta_custom1);
gamma2 = (eta_custom3 - eta_custom2)/(eta_custom3 + eta_custom2);
gamma3 = (eta_air - eta_custom3)/(eta_air + eta_custom3);

d = wavelength10;

gamma2t = (gamma2 + gamma3.*exp(-1j.*2.*beta_custom3.*d))./(1 + gamma2.*gamma3.*exp(-1j.*2.*beta_custom3.*d));
gamma1t = (gamma1 + gamma2t.*exp(-1j.*2.*beta_custom2.*d))./(1 + gamma1.*gamma2t.*exp(-1j.*2.*beta_custom2.*d));
gamma = (gamma0 + gamma1t.*exp(-1j.*2.*beta_custom1.*d))./(1 + gamma0.*gamma1t.*exp(-1j.*2.*beta_custom1.*d));

tau = sqrt(1 - abs(gamma).^2);

% Import HFSS results

hfss9 = readmatrix("S Parameter Plot - dB S11 and S12 Multilayer (2-3-4) Slab.csv");

figure(9)
plot(hz, mag2db(abs(gamma)), hz, mag2db(abs(tau)), hz, hfss9(:,2), hz, hfss9(:,3))
xlabel('Frequency in Hz')
ylabel('dB')
legend('S11','S21','S11 HFSS','S21 HFSS')
title('S Paramater Plot - dB S11 and S21 Multilayer (2-3-4) Slab')

% === 2-8-2 ===
gamma0 = (eta_custom1 - eta_air)/(eta_custom1 + eta_air);
gamma1 = (eta_custom7 - eta_custom1)/(eta_custom7 + eta_custom1);
gamma2 = (eta_custom1 - eta_custom7)/(eta_custom1 + eta_custom7);
gamma3 = (eta_air - eta_custom1)/(eta_air + eta_custom1);

d = wavelength10;

gamma2t = (gamma2 + gamma3.*exp(-1j.*2.*beta_custom1.*d))./(1 + gamma2.*gamma3.*exp(-1j.*2.*beta_custom1.*d));
gamma1t = (gamma1 + gamma2t.*exp(-1j.*2.*beta_custom7.*d))./(1 + gamma1.*gamma2t.*exp(-1j.*2.*beta_custom7.*d));
gamma = (gamma0 + gamma1t.*exp(-1j.*2.*beta_custom1.*d))./(1 + gamma1.*gamma2t.*exp(-1j.*2.*beta_custom1.*d));

tau = sqrt(1 - abs(gamma).^2);

% Import HFSS results

hfss10 = readmatrix("S Parameter Plot - dB S11 and S12 Multilayer (2-8-2) Slab.csv");

figure(10)
plot(hz, mag2db(abs(gamma)), hz, mag2db(abs(tau)), hz, hfss10(:,2), hz, hfss10(:,3))
xlabel('Frequency in Hz')
ylabel('dB')
legend('S11','S21','S11 HFSS','S21 HFSS')
title('S Parameter Plot - dB S11 and S21 Multilayer (2-8-2) Slab')

% === 8-2-8 ===
gamma0 = (eta_custom7 - eta_air)/(eta_custom7 + eta_air);
gamma1 = (eta_custom1 - eta_custom7)/(eta_custom1 + eta_custom7);
gamma2 = (eta_custom7 - eta_custom1)/(eta_custom7 + eta_custom1);
gamma3 = (eta_air - eta_custom7)/(eta_air + eta_custom7);

d = wavelength10;

gamma2t = (gamma2 + gamma3.*exp(-1j.*2.*beta_custom7.*d))./(1 + gamma2.*gamma3.*exp(-1j.*2.*beta_custom7.*d));
gamma1t = (gamma1 + gamma2t.*exp(-1j.*2.*beta_custom1.*d))./(1 + gamma1.*gamma2t.*exp(-1j.*2.*beta_custom1.*d));
gamma = (gamma0 + gamma1t.*exp(-1j.*2.*beta_custom7.*d))./(1 + gamma0.*gamma1t.*exp(-1j.*2.*beta_custom7.*d));

tau = sqrt(1 - abs(gamma).^2);

% Import and plot HFSS results
hfss11 = readmatrix('S Parameter Plot - dB S11 and S12 Multilayer (8-2-8) Slab.csv');

figure(11)
plot(hz, 20.*log10(gamma), hz, 20.*log10(tau), hz, hfss11(:,2), hz, hfss11(:,3))
xlabel('Frequency in Hz')
ylabel('dB')
legend('S11','S21','S11 HFSS','S21 HFSS')
title('S Parameter Plot - dB S11 and S21 Multilayer (8-2-8) Slab')

% === Air - 5 Dielectrics ===

% Apply Theory of Small Reflections

beta_custom1 = omega.*sqrt(mu.*epsilon_custom1);
beta_custom4 = omega.*sqrt(mu.*epsilon_custom4);
beta_custom5 = omega.*sqrt(mu.*epsilon_custom5);
beta_custom6 = omega.*sqrt(mu.*epsilon_custom6);

gamma0 = (eta_custom1 - eta_air)./(eta_custom1 + eta_air);
gamma1 = (eta_custom4 - eta_custom1)./(eta_custom4 + eta_custom1);
gamma2 = (eta_custom5 - eta_custom4)./(eta_custom5 + eta_custom4);
gamma3 = (eta_custom6 - eta_custom5)./(eta_custom6 + eta_custom5);
gamma4 = (eta_custom7 - eta_custom6)./(eta_custom7 + eta_custom6);

d = 0.013;

gamma = gamma0 + gamma1.*exp(-1j.*2.*beta_custom1.*d) + gamma2.*exp(-1j.*2.*(beta_custom1.*d + beta_custom4.*d)) + gamma3.*exp(-1j.*2.*(beta_custom1.*d + beta_custom4.*d + beta_custom5.*d)) + gamma4.*exp(-1j.*2.*(beta_custom1.*d + beta_custom4.*d + beta_custom5.*d + beta_custom6.*d));
tau = sqrt(1 - abs(gamma).^2);

% Import HFSS results

hfss12 = readmatrix("S Parameter Plot - dB S11 and S12 Multilayer (2-3.5-5-6.csv");

figure(12)
plot(hz, mag2db(abs(gamma)), hz, mag2db(abs(tau)), hz, hfss12(:,2), hz, hfss12(:,3))
xlabel('Frequency in Hz')
ylabel('dB')
legend('S11','S21','S11 HFSS','S21 HFSS')
title('S Paramter Plot - dB S11 and S21 Multilayer (2-3.5-5-6.5-8)')



