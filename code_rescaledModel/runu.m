clear all; close all;
% Input parameters.
Lambda = 16;
Nw = 8;
MC = 1e3;
theta = 10;
gibd = 0.24;
fi = 1;
% time step parameters
dt = 5E-4;
nout = 1;
tfin = dt;

% Run the simulation.
J = 2*Lambda;
SolverKdV_SymplecticM4a_MC(J, tfin, dt, nout, MC, 0, Nw, theta, gibd, fi)

% Load the output file.
load('output.mat')
J = size(uarray,1);
MC = size(uarray,2);
Ntsteps = size(uarray,3);
% Flatten the array of u into a single list.
ulist = reshape(uarray,1,[]);
dulist = reshape(duarray,1,[]);
% Write output to a text file.
fileID = fopen('ulist.txt','w');
fprintf(fileID,'# Number of grid points, trajectories, and time steps.\n');
fprintf(fileID,'%d\n',[J,MC,Ntsteps]);
fprintf(fileID,'# Values of u\n');
fprintf(fileID,'%12.8f\n', ulist);
fprintf(fileID,'# Values of du\n');
fprintf(fileID,'%12.8f\n', dulist);
fclose(fileID);

% Make a plot to see the waves.
%x=(-pi:2*pi/J:pi-2*pi/J)';
%contourf(T,x,uarray(:,1,:));

% Make a histogram of u.
figure(1); histogram(ulist);
figure(2); histogram(ulist); set(gca,'yscale','log');
skewu = skewness(ulist)