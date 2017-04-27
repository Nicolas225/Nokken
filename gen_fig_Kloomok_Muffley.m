clear; close all; clc


%% input
R0 = 75; %pitch radius in mm
start_lift = 30; %start lift in mm
end_lift = 0; %end lift in mm
cycloidtype = 6; %type of cycloide

%% calculations
beta_vec = 5:200;
rho_min = 0*beta_vec;

L0 = start_lift;
L1 = end_lift;
for i = 1:length(beta_vec)
  beta = beta_vec(i)*pi/180;
  theta = linspace(0, beta_vec(i), 100)*pi/180;
  
  if cycloidtype==1
    L=L1-L0;
    S=L*(theta/beta-sin(pi*theta/beta)/pi)+L0;
    V=L/beta*(1-cos(pi*theta/beta));
    A=pi*L/beta^2*(sin(pi*theta/beta));
    
  elseif cycloidtype==2
    L=L1-L0;
    S=L*(theta/beta+sin(pi*theta/beta)/pi)+L0;
    V=L/beta*(1+cos(pi*theta/beta));
    A=-pi*L/beta^2*(sin(pi*theta/beta));
    
  elseif cycloidtype==3
    L=-(L1-L0);
    S=L*(1-theta/beta+sin(pi*theta/beta)/pi)+L0-L;
    V=-L/beta*(1-cos(pi*theta/beta));
    A=-pi*L/beta^2*(sin(pi*theta/beta));
    
  elseif cycloidtype==4
    L=-(L1-L0);
    S=L*(1-theta/beta-sin(pi*theta/beta)/pi)+L0-L;
    V=-L/beta*(1+cos(pi*theta/beta));
    A=pi*L/beta^2*(sin(pi*theta/beta));
    
  elseif cycloidtype==5
    L=L1-L0;
    S=L*(theta/beta-sin(2*pi*theta/beta)/2/pi)+L0;
    V=L/beta*(1-cos(2*pi*theta/beta));
    A=2*pi*L/beta^2*(sin(2*pi*theta/beta));
    
  elseif cycloidtype==6
    L=-(L1-L0);
    S=L*(1-theta/beta+sin(2*pi*theta/beta)/2/pi)+L0-L;
    V=-L/beta*(1-cos(2*pi*theta/beta));
    A=-2*pi*L/beta^2*(sin(2*pi*theta/beta));
  else
    error('cycloid type must be a number between 1 and 6')
  end
  
  rho = ((R0+S).^2 + V.^2).^(3/2) ./ ((R0+S).^2 + 2*V.^2 - (R0+S).*A);
  rho_min(i) = min(rho(rho>0));
  
end

figure;
plot(beta_vec, rho_min);
grid
xlabel('beta (degrees)')
ylabel('rho min (mm)')
