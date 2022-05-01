clc
clear all 
format long g

%% Main Program
clc
clear all
format long g
% Input Initial value
iframe = [
1400.966 -1029.940;
6488.608 -7930.481;
3792.613 -8116.757;
9062.577 -8805.310;
%5988.875 -3726.625;
1111.500 -9706.804
];

mframe = [
558667.793 1639796.997;
562246.621 1642483.679;
562361.250 1641049.735;
562696.450 1643849.215;
%560030.781 1642197.631;
563213.131 1639649.640
];
%%
clc
% Affine Transformation 
% A matrix function at bottom of this file
A = AffineMatrixACreation(iframe,mframe)
% L matrix
L  = AffineLmatrix(mframe)
% N and U matrix
N = A'*A
U = A'*L
% X
X = inv(N)*U
%RMSE
m = length(L);
n = 6; % for affine transformation
V = A*X-L;
RMSE = sqrt((V'*V)/(m-n))
% Visualize Value
coordinate = [5988.875 -3726.625];
calculateCoordinate = affineENCoordinate(iframe,mframe,coordinate)
referenceCoordinate = [560030.781 1642197.631]
dE = calculateCoordinate(1)-referenceCoordinate(1)
dN = calculateCoordinate(2)-referenceCoordinate(2)
dS = sqrt(dE^2+dN^2)
RMSE
diff = [dE dN]
%%
clc
% Helmert Transformation
% A matrix function at bottom of this file
A = HelmertMatrixACreation(iframe,mframe)
% L matrix
L = HelmertLmatrix(mframe)
% X
X = HelmertTransformation(iframe,mframe)
a = X(1);
b = X(2);
theta = rad2deg(atan(b/a))
m_scale = a/cosd(theta)
%RMSE
m = length(L);
n = 4; % for helmert transformation
V = A*X-L;
RMSE = sqrt((V'*V)/(m-n))
coordinate = [5988.875 -3726.625];
calculateCoordinate2 = helmertENCoordinate(iframe,mframe,coordinate)
referenceCoordinate2 = [560030.781 1642197.631]
dE2 = calculateCoordinate2(1)-referenceCoordinate2(1)
dN2 = calculateCoordinate2(2)-referenceCoordinate2(2)
dS2 = sqrt(dE2^2+dN2^2)
theta
m_scale
RMSE
diff = [dE2 dN2]

%% Function

% AFFINE TRANSFORMATION
% this function is to Create an A matrix for Affine transformation
function A = AffineMatrixACreation(iframe,mframe)
    A = [];
    for i = 1 : length(iframe)
        each = [
    iframe(i,1) iframe(i,2) 1 0 0 0;
    0 0 0 iframe(i,1) iframe(i,2) 1;
    ];
        A = [A;each];
    end
end

% this function is to Create an L matrix for Affine transformation
function L = AffineLmatrix(mframe)
    L = [];
    for i = 1 : length(mframe)
        for j = 1:2
            val = mframe(i,j);
            L = [L;val];
        end
    end
end
function X = AffineTransformation(iframe,mframe)
    A = AffineMatrixACreation(iframe,mframe);
    % L matrix
    L = AffineLmatrix(mframe);
    % N and U matrix
    N = A'*A;
    U = A'*L;
    % X
    X = inv(N)*U;
end

% this function is a function to calculate EN coordinate (mframe
% coordinate) from psuedo parameter that we calculate by leastsquared
function EN = affineENCoordinate(iframe,mframe,coordinate)
    X = AffineTransformation(iframe,mframe);
    a = X(1);
    b = X(2);
    c = X(3);
    d = X(4);
    e = X(5);
    f = X(6);
    mx = coordinate(1);
    my = coordinate(2);
    E = mx*a + my*b + c;
    N = mx*d + my*e + f;
    EN = [E N];
end

% HELMERT TRANSFORMATION
function A = HelmertMatrixACreation(iframe,mframe)
    A = [];
    for i = 1 : length(iframe)
        each = [
    iframe(i,1) iframe(i,2) 1 0 ;
    iframe(i,2) -1*iframe(i,1) 0 1;
    ];
        A = [A;each];
    end
end
% this function is to Create an L matrix for Affine transformation
function L = HelmertLmatrix(mframe)
    L = [];
    for i = 1 : length(mframe)
        for j = 1:2
            val = mframe(i,j);
            L = [L;val];
        end
    end
end
function X = HelmertTransformation(iframe,mframe)
    A = HelmertMatrixACreation(iframe,mframe);
    % L matrix
    L = HelmertLmatrix(mframe);
    % N and U matrix
    N = A'*A;
    U = A'*L;
    % X
    X = inv(N)*U;
end
function EN = helmertENCoordinate(iframe,mframe,coordinate)
    X = HelmertTransformation(iframe,mframe);
    a = X(1)
    b = X(2)
    c = X(3)
    d = X(4)
    mx = coordinate(1);
    my = coordinate(2);
    E = mx*a + my*b + c;
    N = mx*(-1)*b + my*a + d;
    EN = [E N];
end