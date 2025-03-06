function [SYS,dSYS] = generateSateSpace(Par,Ts)
% Function that generates the state space model (continuous and discrete
% time) of a one storage building with an AMD on the top floor
%
% Inputs: 
% - Parameters = [k_f,c_f,m_f,m_c,c_c]
%       Parameter vector with:
%           k_f -> spring coefficient last floor
%           c_f -> damper coefficient last floor
%           m_f -> mass last floor
%           m_c -> mass car (AMD)
%           c_c -> damper (friction) car-building
% - Ts: Sample time used in the discretization
%
% Outputs
% SYS -> continuous time state space
% dSYS -> discrete time state space model using a zero order holder 
%
% State-Space Matrices
k_f = Par(1);
c_f = Par(2);
m_f = Par(3);
m_c = Par(4);
c_c = Par(5);

A = [ 0                                        1     0                   0;
     -k_f/m_f                           -c_f/m_f     0             c_c/m_f;
      0                                        0     0                   1;
      k_f/m_f                            c_f/m_f     0     -c_c/m_c-c_c/m_f];

B = [ 0    0;
     -1   0;
      0    0;
      0   1/m_c];

C = [-k_f/m_f -c_f/m_f  0  c_c/m_f;
            0        0  1        0];

D = [0   0;
    0   0]; % [-1   0];

SYS = ss(A,B,C,D);
dSYS = c2d(SYS,Ts,'zoh');
end