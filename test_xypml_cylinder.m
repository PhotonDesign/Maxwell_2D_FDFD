%The computational region is a square of size xlength*ylength, i.e. Dim(1)*Dim(2), with pml of
%'thickness' enclosed both in x and y.

%Any function (source, epsilon, solution) is defined on the M*N grid points

%The ordering of the grid goes like [1 2 ... N; N+1 N+2 ... 2N; ...] 


%close all;

Dim = [2 2];
h = 0.01;

N = round(Dim(1)/h);%num of x dim grid points
M = round(Dim(2)/h);%num of y dim grid points
omega = 2*pi/0.5;
thickness = 0.2; %pml thickness
beta = 10; %pml strength
BC = {{'pml', [thickness,beta]}, {'pml', [thickness,beta]}};

Matx = ones(M,N); %mux, staggered
Maty = ones(M,N); %muy, staggered
const_bg = 1;
Matz = const_bg*ones(M,N); %epsz, not staggered

R  = 0.3; const = 6;

%define a sphere
for i = 1:M
    for j = 1:N
        if (i - M/2)^2 + (j - N/2)^2 < (R/h)^2
            Matz(i,j) = const;
        end
    end
end

%imagesc(Matz);


%source
Source = zeros(M,N);
Source(round(M/2),round(N/4)) = -1i*omega/h^2;

Solution = Scattering_Solve(omega,Dim,h,BC,Matx,Maty,Matz, Source);

figure;
imagesc([0 Dim(2)], [0 Dim(1)],Matz);



%
figure;
imagesc([0 Dim(2)], [0 Dim(1)], real(Solution));
figure;
imagesc([0 Dim(2)], [0 Dim(1)],imag(Solution));



%real part check with analytical solution

x  = linspace(0,Dim(1),N);
y  = linspace(0,Dim(2),M);
%Matz = const*ones(M,N); %epsz, not staggered

numer_alpham = @(m,R,epsilon,k) besselj(m,k*R)*(besselj(m-1,sqrt(epsilon)*k*R) - besselj(m+1,sqrt(epsilon)*k*R))*0.5 ...
    - (1/sqrt(epsilon))*(besselj(m-1,k*R) - besselj(m+1,k*R))*0.5*besselj(m,sqrt(epsilon)*k*R);
denorm_alpham = @(m,R,epsilon,k) besselh(m,k*R)*(besselj(m-1,sqrt(epsilon)*k*R) - besselj(m+1,sqrt(epsilon)*k*R))*0.5 ...
    - (1/sqrt(epsilon))*(besselh(m-1,k*R) - besselh(m+1,k*R))*0.5*besselj(m,sqrt(epsilon)*k*R);


numer_betam = @(m,R,epsilon,k) (1/sqrt(epsilon))*besselj(m,k*R)*(besselh(m-1,k*R) - besselh(m+1,k*R))*0.5 ...
    - (1/sqrt(epsilon))*(besselj(m-1,k*R) - besselj(m+1,k*R))*0.5*besselh(m,k*R);


alpham = @(m,R,epsilon,k) -numer_alpham(m,R,epsilon,k)/denorm_alpham(m,R,epsilon,k);
betam = @(m,R,epsilon,k) -numer_betam(m,R,epsilon,k)/denorm_alpham(m,R,epsilon,k);



add_thm = zeros(1,N);
for  i = 1:N
r = abs(h*(i - N/2)); rp =  h*(N/2 - N/4);
t = pi; tp = pi;


m = 20;
for  j = -m:m
    if i < N/4
        add_thm(i) = add_thm(i) + besselh(j,sqrt(const_bg)*omega*r)*besselj(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(t-tp))...
            + alpham(j,R,const,omega)*besselh(j,sqrt(const_bg)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(t-tp));
    elseif i >= N/4 && i < (N/2 - R/h)
        add_thm(i) = add_thm(i) + besselj(j,sqrt(const_bg)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(t-tp))...
        + alpham(j,R,const,omega)*besselh(j,sqrt(const_bg)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(t-tp));
    elseif i < (N/2) && i >= (N/2 - R/h)
        %add_thm(i) = 0 ;
        add_thm(i) = add_thm(i) + ... %besselj(j,sqrt(const_bg)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(t-tp))...
        + betam(j,R,const,omega)*besselj(j,sqrt(const)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(t-tp));
    elseif i >= (N/2) && i <= (N/2 + R/h)
        %add_thm(i) = 0 ;
        add_thm(i) = add_thm(i) + ... %besselj(j,sqrt(const_bg)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(-tp))...
        + betam(j,R,const,omega)*besselj(j,sqrt(const)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(-tp));
    elseif i > (N/2 + R/h) && i < (3*N/4)
        add_thm(i) = add_thm(i) + besselj(j,sqrt(const_bg)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(-tp))...
        + alpham(j,R,const,omega)*besselh(j,sqrt(const_bg)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(-tp));
    else
        add_thm(i) = add_thm(i)+besselh(j,sqrt(const_bg)*omega*r)*besselj(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(-tp))...
        + alpham(j,R,const,omega)*besselh(j,sqrt(const_bg)*omega*r)*besselh(j,sqrt(const_bg)*omega*rp)*exp(1i*j*(-tp));
    end
    
%     if r == rp
%         add_thm(i) = 0;
%     end
    
end
% disp(add_thm(i));
% result2 = besselh(0,sqrt(const)*omega*norm([cos(t)*r;sin(t)*r] - [cos(tp)*rp;sin(tp)*rp]));
% disp(result2);
% disp(' ');
end

% figure;
% plot(imag(add_thm));


figure;
plot(x,imag(Solution(round(M/2),:)))
hold on;
x  = linspace(0,Dim(1),N);
plot(x,imag(1i*(-1i)*omega*add_thm/4),'--');
hold off;
legend('FD solution','addition formula');
ylim([-10 10]);

figure;
plot(x,real(Solution(round(M/2),:)))
hold on;
x  = linspace(0,Dim(1),N);
plot(x,real(1i*(-1i)*omega*add_thm/4),'--');
hold off;
legend('FD solution','addition formula');
ylim([-10 10]);

%%

