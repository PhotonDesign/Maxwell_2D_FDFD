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
Matz = const*ones(M,N); %epsz vacuum, not staggered

Source = zeros(M,N);
Source(round(M*0.7),:) = 1/h; %a line source

Vac_Sol = Scattering_Solve(omega,Dim,h,BC,Matx,Maty,Matz, Source); %vacuum solution


%define the slab
const = 12; %epsz vacuum slab refraction index
n = 50;
Layer_Thickness = linspace(0.1,0.3,n); %thickness of the slab to be varied

%plane position to measure flux
T_plane = 0.2;
R_plane = 0.8;


%R and T calculation
R = zeros(1,n);
T = zeros(1,n);
for i= 1:n
    Matz(round(M*0.3):round(M*(Layer_Thickness(i)+0.3)),:) = const;
    Sol = Scattering_Solve(omega,Dim,h,BC,Matx,Maty,Matz, Source);
    difference = Sol(:,round(N/2))-Vac_Sol(:,round(N/2));
    T(i) = abs(Sol(round(T_plane*M),round(N/2))/Vac_Sol(round(T_plane*M),round(N/2)))^2;
    R(i) = abs(difference(round(R_plane*M))/Vac_Sol(round(R_plane*M),round(N/2)))^2;
end

%analytical solution
r = -(1-sqrt(const))/(1+sqrt(const));
phi = omega*sqrt(const)*(Layer_Thickness*Dim(2)+h);
R_exact = abs((r-r*exp(-2*1i*phi))./(1-r^2*exp(-2*1i*phi))).^2;

%plot T and R
plot(Layer_Thickness,R);
hold on;
plot(Layer_Thickness,T);
plot(Layer_Thickness,T+R);
plot(Layer_Thickness,R_exact);
plot(Layer_Thickness,1-R_exact);
legend('R','T','R+T', 'R exact', 'T exact');
xlabel('layer thickness');