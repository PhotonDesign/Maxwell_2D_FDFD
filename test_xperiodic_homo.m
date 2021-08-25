%be careful about reshape [N M] and Maxtrix(:) operation
Dim = [4 4];
h = 0.01;

N = round(Dim(1)/h);%num of x dim grid points
M = round(Dim(2)/h);%num of y dim grid points
omega = 2*pi/0.5;
thickness = 0.6;
beta = 10;
BC = {{'periodic'}, {'pml', [thickness,beta]}};

Matx = ones(M,N); %mux, staggered
Maty = ones(M,N); %muy, staggered
const = 1;
Matz = const*ones(M,N); %epsz, not staggered

Source = zeros(M,N);
Source(round(M/2),:) = -1i*omega/h;

Solution = Scattering_Solve(omega,Dim,h,BC,Matx,Maty,Matz, Source);

figure;
imagesc(real(Solution));
figure;
imagesc(imag(Solution));


%real part check
figure;
plot(1:M,real(Solution(:,round(N/2))));
hold on;
y = linspace(0,Dim(2),M);
plot(1:M,omega*cos(omega*sqrt(const)*abs(y-y(round(M/2))))*1/(2*omega*sqrt(const)));
hold off;
legend('FD solution', 'exact solution');

%imag part check
figure;
plot(1:M,imag(Solution(:,round(N/2))));
hold on;
y = linspace(0,Dim(2),M);
plot(1:M,omega*sin(omega*sqrt(const)*abs(y-y(round(M/2))))*1/(2*omega*sqrt(const)));
hold off;
legend('FD solution', 'exact solution');