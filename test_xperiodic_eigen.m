%The computational region is a square of size xlength*ylength, i.e. Dim(1)*Dim(2), with pml of
%'thickness' enclosed both in x and y.

%Any function (source, epsilon, solution) is defined on the M*N grid points

%The ordering of the grid goes like [1 2 ... N; N+1 N+2 ... 2N; ...] 

%This is an eigenproblem.
%Eigen_Maxwell is the eps^(-1)curl mu^(-1) curl operator.
%The user can obtain eigenmodes by eigensolvers in Matlab

%pml size and strength should be chosen carefully to unveil QNMs and the
%continuum in a distinctive manner.



Dim = [1 12];
h = 0.01;

N = round(Dim(1)/h);%num of x dim grid points
M = round(Dim(2)/h);%num of y dim grid points
thickness = 0.25*Dim(2); %pml thickness
beta = 20; %pml strength
BC = {{'periodic'}, {'pml', [thickness,beta]}};

Matx = ones(M,N); %mux, staggered
Maty = ones(M,N); %muy, staggered
Matz = ones(M,N); %epsz, not staggered
const = 12;
d = 2;
Matz(round(0.4*Dim(2)/h):round((0.4*Dim(2)+d)/h),:) = const;

%Eigen_Maxwell is the eps^(-1)curl mu^(-1) curl operator.
Eigen_Maxwell = Eigen_Operator(Dim,h,BC,Matx,Maty,Matz);
%D: w^2, V: mode profile, flag: 0 if converge
[V,D,flag] = eigs(Eigen_Maxwell, 60, 'smallestabs'); %flag = 0 eigenvalue convergences


D = diag(sqrt(D));
figure;
plot(real(D), imag(D), '*');
hold on;
%analytical
w = [(pi + 1i*log((sqrt(12)-1)/(sqrt(12)+1))), (2*pi + 1i*log((sqrt(const)-1)/(sqrt(const)+1)))]/(d*sqrt(const));
plot(real(w), imag(w), 's');
hold off;
legend('FDFD solver', 'analytic');
xlim([0 2])

%find(D>0.2)
eig_vec = reshape(V(:,28),[N,M]).';
figure;
imagesc(real(eig_vec));
figure;
plot(real(eig_vec(:,round(N/2))));

eig_vec = reshape(V(:,55),[N,M]).';
figure;
imagesc(real(eig_vec));
figure;
plot(real(eig_vec(:,round(N/2))));

