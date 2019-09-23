clear all; close all;

% Model parameters from experiments
epsi0 = 0.017;  % amplitude-to-depth ratio
del0 = 0.22;    % depth-to-wavelength ratio
Drat = 0.24;    % depth ratio.
% Simulation parameters
Lambda = 16;    % Reference 16
Nw = Lambda/2;  % Reference Lambda/2
MC = 1E2;       % Reference 1E4
% time step parameters
dt = 5E-4;  % Reference 5E-4
nout = 100; % Reference 100
tfin = 5.; % Reference 10

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

% Run the simulation.
tic;
C2 = (2/3)*pi^2*del0/Nw^2
C3 = 3/2*pi^(1/2)*epsi0/del0
SolverKdV_SymplecticM4a_MC(C2,C3,Drat,Lambda,MC,theta,gibd,fi,dt,nout,tfin);
cputime = round(toc/60, 2, 'significant')

% Load the output file.
load('output.mat')
JJ = size(uarray,1);
MC = size(uarray,2);
Ntsteps = size(uarray,3);
% Flatten the array of u into a single list.
ulist = reshape(uarray,1,[]);
dulist = reshape(duarray,1,[]);
% Write output to a text file.
fileID = fopen('ulist.txt','w');
% Input parameters
fprintf(fileID,'#Basic integers: JJ, MC, Nt, Nw.\n');
out1 = [JJ,MC,nout,Nw]; fprintf(fileID,'%d\n', out1);
fprintf(fileID,'#Basic floats: theta, tfin, cputime (mins).\n');
out2 = [theta, tfin, cputime]; fprintf(fileID,'%9.5f\n', out2);
fprintf(fileID,'# All other inputs: C2,C3,Drat,gibd,fi,dt.\n');
out3 = [C2,C3,Drat,gibd,fi,dt]; fprintf(fileID,'%9.5f\n', out3);
fprintf(fileID,'# Pad with zeros until index 20.\n');
len = 20 - length(out1)-length(out2)-length(out3);
fprintf(fileID,'%d\n', zeros(len,1) );
% Microstates
fprintf(fileID,'# Values of u\n'); fprintf(fileID,'%9.5f\n', ulist);
fprintf(fileID,'# Values of du\n'); fprintf(fileID,'%9.5f\n', dulist);
fclose(fileID);

% Make a plot to see the waves.
%x=(-pi:2*pi/J:pi-2*pi/J)';
%contourf(T,x,uarray(:,1,:));

% Make a histogram of u.
figure(1); histogram(ulist);
figure(2); histogram(ulist); set(gca,'yscale','log');
skewu = skewness(ulist)