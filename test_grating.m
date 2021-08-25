ulen = 350; %unit length as 350 nm
Dim = [1, 8]; %Dim = [x dimension, y dimension]
h = 0.01; %grid size
N = round(Dim(1)/h);%num of x dim grid points
M = round(Dim(2)/h);%num of y dim grid points
thickness = 0.6; beta = 10; %PML parameter
BC = {{'periodic'},{'pml', [thickness,beta]}};

eps_GaN = 2.47^2;
Matx = ones(M,N); %mux, staggered
Maty = ones(M,N); %muy, staggered
Matz = eps_GaN*ones(M,N); %epsz vacuum, not staggered

Source = zeros(M,N);
Source(round(M*0.2),:) = 1/h; %a line source to create plane waves

lambda = 480/ulen; w = 2*pi/lambda; %c = 1
%% Vacuum solution
Ez_vac = Scattering_Solve(w,Dim,h,BC,Matx,Maty,Matz, Source); %vacuum solution
%get H field
[Hx_vac,~] = Get_Comp_Fields(w, BC, Dim,h, Ez_vac);
flux_plane = 0.3*Dim(2); %measure power flow at y=flux_plane
Poynting_vac = 0.5*real(conj(Hx_vac).*Ez_vac); %vacuum poynting vector
P_vac = trapz(Poynting_vac(round(flux_plane/h),:))*h; %total incoming power in vacuum
%% Simulation setup, creating a block
eps_vac = 1; %epsz vacuum slab refraction index
Layer_Thickness = 300/ulen; %thickness of the slab to be varied
Block_Size = 0.5*350/ulen; %size of the block
%grid position of the block
block_posy = round(0.5*Dim(2)/h):round((0.5*Dim(2)+Layer_Thickness)/h); block_posx = 1:round(Block_Size/h);
Matz(block_posy,block_posx) = eps_vac;
%check
imagesc(Matz);
%% Solve for the problem
Ez = Scattering_Solve(w,Dim,h,BC,Matx,Maty,Matz, Source); 
Ez_scat = Ez - Ez_vac; %get the scattering field only to calculate reflection
%get H field 
[Hx_scat,~] = Get_Comp_Fields(w, BC, Dim,h, Ez_scat);
%% Computing power flow
Poynting = 0.5*real(conj(Hx_scat).*Ez_scat);
R = trapz(Poynting(round(flux_plane/h),:))*h;
reflection = -R/P_vac;
disp('the relfection is:');
disp(reflection);
%% Visualizing fields
figure; imagesc(real(Ez)); figure; plot(real(Ez(:,10))); %check total field
figure; imagesc(real(Ez_scat)); figure; plot(real(Ez_scat(:,10))); %check scattered field