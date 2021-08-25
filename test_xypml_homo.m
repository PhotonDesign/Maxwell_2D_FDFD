%The computational region is a square of size xlength*ylength, i.e. Dim(1)*Dim(2), with pml of
%'thickness' enclosed both in x and y.

%Any function (source, epsilon, solution) is defined on the M*N grid points

%The ordering of the grid goes like [1 2 ... N; N+1 N+2 ... 2N; ...] 

Dim = [2 2];
h = 0.01;

N = round(Dim(1)/h);%num of x dim grid points
M = round(Dim(2)/h);%num of y dim grid points
omega = 2*pi/0.5;
thickness = 0.4; %pml thickness
beta = 10; %pml strength
BC = {{'pml', [thickness,beta]}, {'pml', [thickness,beta]}};

Matx = ones(M,N); %mux, staggered
Maty = ones(M,N); %muy, staggered
const = 1;
Matz = const*ones(M,N); %epsz, not staggered

%source
Source = zeros(M,N);
Source(round(M/2),round(N/2)) = -1i*omega/h^2;

Solution = Scattering_Solve(omega,Dim,h,BC,Matx,Maty,Matz, Source);

figure;
imagesc(real(Solution));
figure;
imagesc(imag(Solution));

%real part check with analytical solution

%x direction
figure;
plot(1:N,real(Solution(round(M/2),:)))
hold on;
x  = linspace(0,Dim(1),N);
plot(1:N,real(1i*(-1i)*omega*besselh(0,1,sqrt(const)*omega*abs(x-x(round(N/2))))/4));
hold off;
legend('FD solution', 'exact solution');
%y direction
figure;
plot(1:M,real(Solution(:,round(N/2))))
hold on;
y  = linspace(0,Dim(2),M);
plot(1:M,real(1i*(-1i)*omega*besselh(0,1,sqrt(const)*omega*abs(y-y(round(M/2))))/4));
hold off;
legend('FD solution', 'exact solution');

%imag part check with analytical solution

%x direction
figure;
plot(1:N,imag(Solution(round(M/2),:)))
hold on;
x  = linspace(0,Dim(1),N);
plot(1:N,imag(1i*(-1i)*omega*besselh(0,1,sqrt(const)*omega*abs(x-x(round(N/2))))/4));
hold off;
legend('FD solution', 'exact solution');
%y direction
figure;
plot(1:M,imag(Solution(:,round(N/2))))
hold on;
y  = linspace(0,Dim(2),M);
plot(1:M,imag(1i*(-1i)*omega*besselh(0,1,sqrt(const)*omega*abs(y-y(round(M/2))))/4));
hold off;
legend('FD solution', 'exact solution');


