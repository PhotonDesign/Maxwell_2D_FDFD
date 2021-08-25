%The computational region is a square of size xlength*ylength, i.e. Dim(1)*Dim(2), with pml of
%'thickness' enclosed both in x and y.

%Any function (source, epsilon, solution) is defined on the M*N grid points

%The ordering of the grid goes like [1 2 ... N; N+1 N+2 ... 2N; ...] 

%This is an eigenproblem.
%Eigen_Maxwell is the eps^(-1)curl mu^(-1) curl operator.
%The user can obtain eigenmodes by eigensolvers in Matlab

%pml size and strength should be chosen carefully to unveil QNMs and the
%continuum in a distinctive manner.



Dim = [4 4];
h = 0.01;

N = round(Dim(1)/h);%num of x dim grid points
M = round(Dim(2)/h);%num of y dim grid points
thickness = 1; %pml thickness
beta = 20; %pml strength
BC = {{'pml', [thickness,beta]}, {'pml', [thickness,beta]}};

Matx = ones(M,N); %mux, staggered
Maty = ones(M,N); %muy, staggered
Matz = ones(M,N); %epsz, not staggered
const = 12;
Matz(round(0.4*Dim(1)/h):round(0.6*Dim(1)/h),round(0.4*Dim(2)/h):round(0.6*Dim(2)/h)) = const;

%Eigen_Maxwell is the eps^(-1)curl mu^(-1) curl operator.
Eigen_Maxwell = Eigen_Operator(Dim,h,BC,Matx,Maty,Matz);
%D: w^2, V: mode profile, flag: 0 if converge
[V,D,flag] = eigs(Eigen_Maxwell, 200, 'smallestabs'); %flag = 0 eigenvalue convergences


D = diag(sqrt(D));
figure;
plot(real(D), imag(D), '*');

eig_vec = reshape(V(:,165),[N,M]).';
figure;
imagesc(real(eig_vec));
figure;
plot(real(eig_vec(:,round(N/2))));

