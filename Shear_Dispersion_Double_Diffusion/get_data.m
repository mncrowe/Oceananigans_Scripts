clear
close all

filename = "equal_diff";

DS = 1;
Dt = 2;

S = ncread(filename + "_tracers.nc", "S");
T = ncread(filename + "_tracers.nc", "T");

u = ncread(filename + "_velocity.nc", "u");
w = ncread(filename + "_velocity.nc", "w");

t = ncread(filename + "_velocity.nc", "time");
x = ncread(filename + "_velocity.nc", "x_caa");
z = ncread(filename + "_velocity.nc", "z_aac");

S = S + x * DS / x(end);
T = T + x * Dt / x(end);

B = T - S;

splot(S(:,:,1), x, z, xlabel='x', ylabel='z', title='S(t = 0)')
splot(T(:,:,1), x, z, xlabel='x', ylabel='z', title='T(t = 0)')
splot(B(:,:,1), x, z, xlabel='x', ylabel='z', title='B(t = 0)')

splot(squeeze(B(:, 51, :)), x, t, xlabel='x', ylabel='t', title='B(x, 0, t)')