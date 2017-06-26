clear all; close all;
%% setup the parameters
df = 0.01;
v=3;
omega = [0:df:20]*2*pi;
lambda = 2*pi*v/max(omega);
Q=120;
rmin = 200*lambda; rmax = 250*lambda; % distance range for noise source
% x1(0,0), x2(50lambda,0),x3(50lambda,pi) % location of three receivers
N = 10000; % # of noise source
r = rmin + (rmax-rmin).*rand(N,1); % get the random r for noise source
theta = pi*rand(N,1).*sign(rand(N,1)*2-1)+pi*(sign(rand(1)*2-1)>0); % random theta for noise source

figure(1)
polar(theta,r,'.k');
%% calculate the waveform on each receiver in frequency domain

uf1 = zeros(size(omega)); uf2 = uf1; uf3 = uf1;
UF1 = zeros(size([0:df:50-df])); UF2 = UF1; UF3=UF1; % pad zeros to frequency 20-50
for i=1:N
    phi = rand(1,length(omega))*2*pi; % random phase
    r2 = sqrt((r(i)*cos(theta(i))-50*lambda)^2+(r(i)*sin(theta(i)))^2); % distance from x2 to noise source
    r3 = sqrt((r(i)*cos(theta(i))+50*lambda)^2+(r(i)*sin(theta(i)))^2); % distance from x1 to noise source
    uf1 = uf1 + 1/sqrt(r(i))*exp(-1i*omega*r(i)/v).*exp(-omega/v*r(i)/2/Q).*exp(1i*phi);
    uf2 = uf2 + 1/sqrt(r2)*exp(-1i*omega*r2/v).*exp(-omega/v*r2/2/Q).*exp(1i*phi);
    uf3 = uf3 + 1/sqrt(r3)*exp(-1i*omega*r3/v).*exp(-omega/v*r3/2/Q).*exp(1i*phi);
end


UF1(1:length(uf1)) = uf1;
UF2(1:length(uf2)) = uf2;
UF3(1:length(uf3)) = uf3;

% from frequency domain to time domain
U1 = ifft(UF1,'symmetric');
U2 = ifft(UF2,'symmetric');
U3 = ifft(UF3,'symmetric');

Nt = length(U1);
dt = 1/(df*Nt);
T=[0:dt:(Nt-1)*dt];
[b,a] = butter(4,[(2*4*dt),(2*12*dt)],'bandpass');

U1 = filtfilt(b,a,U1);
U2 = filtfilt(b,a,U2);
U3 = filtfilt(b,a,U3);

%% deconvolution interferometry

figure  % plot the waveforms
plot(T,filtfilt(b,a,U1))
hold on
plot(T,filtfilt(b,a,U2))
plot(T,filtfilt(b,a,U3))

UF1 = (fft(U1,length(U1)*2-1));
UF2 = (fft(U2,length(U2)*2-1));
UF3 = (fft(U3,length(U3)*2-1));
G21f = UF1.*conj(UF2)./(abs(UF2).^2+0.0001);
G12f = UF2.*conj(UF1)./(abs(UF1).^2+0.0001);
G32f = UF2.*conj(UF3)./(abs(UF3).^2+0.0001);
G23f = UF3.*conj(UF2)./(abs(UF2).^2+0.0001);
figure
plot([-fliplr(T),T(2:end)],filtfilt(b,a,(real(ifftshift(ifft(G21f))))))
figure
plot([-fliplr(T),T(2:end)],filtfilt(b,a,(real(ifftshift(ifft(G12f))))))