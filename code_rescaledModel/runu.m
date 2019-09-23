clear all; close all;
% Input parameters.
Lambda = 16;
Nw = 8;
MC = 1e3;
% Choose upstream or downstream
down = false;
if down == true
    theta = 13;
    gibd = 0.24;
    fi = 1;
else
    theta = 15;
    gibd = 1.;
    fi = 0;
end
% time step parameters
dt = 5E-4;
nout = 100; % Default 100?
tfin = 10;


% Run the simulation.
J = 2*Lambda; tic;
SolverKdV_SymplecticM4a_MC(J, tfin, dt, nout, MC, 0, Nw, theta, gibd, fi)
cputime = toc

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
fprintf(fileID,'# Number of grid points, trajectories, time steps, ');
fprintf(fileID,'value of Nw, theta, tfin.\n');
fprintf(fileID,'%d\n',[J,MC,Ntsteps,Nw,theta,tfin]);
fprintf(fileID,'# Values of u\n');
fprintf(fileID,'%9.5f\n', ulist);
fprintf(fileID,'# Values of du\n');
fprintf(fileID,'%9.5f\n', dulist);
fclose(fileID);

% Make a plot to see the waves.
%x=(-pi:2*pi/J:pi-2*pi/J)';
%contourf(T,x,uarray(:,1,:));

% Make a histogram of u.
figure(1); histogram(ulist);
figure(2); histogram(ulist); set(gca,'yscale','log');
skewu = skewness(ulist)