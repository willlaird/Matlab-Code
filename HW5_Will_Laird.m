%%
clc
clear all;
%% Initializing Variables

g = 9.81; %% graviational constant (m/s^2).
S = 124; %% wing span m^2.
b = 36; %% wing span in m.
Wmax = 540e3; %% max weight in N.
Wfuel = 180e3; %% weight of fuel in N.
e_prime = 0.8; 
Cd_0 = 0.04; %% CD.

c_t = 1e-4; %%specific fuel consumption (N/Ns).
Ta_sl = 200e3; %%combined thrust from both turbofans in N.
rho_sl = 1.225; %% sea level air density (kg/m^3).
AR = (36^2)/(124);

time = linspace(0,15e3,15e3+1); %%array for time values.

dt = time(2) - time(1); %%time step between points.

%intializing arrays for calculations:
distance = zeros(length(time),1); %%linspace of distances.
Cl = zeros(length(time),1);
Cd = zeros(length(time),1);
rho = zeros(length(time),1); %% density
W = zeros(length(time),1); %%weight
theta = zeros(length(time),1); %%climb angle
fuel_used = zeros(length(time),1);
Tr = zeros(length(time),1); %thrust required.
Ta = zeros(length(time),1); %thrust applied.
Pa = zeros(length(time),1); %%power available.
Pr = zeros(length(time),1); %%power required.
h = zeros(length(time),1); %%altitude
V = zeros(length(time),1); %%Air velocity
V_v = zeros(length(time),1); %%veritcal velocity
V_h = zeros(length(time),1); %% horizontal velocity



%% Initial conditions for the climb:

W(1) = Wmax;

Ta(1) = Ta_sl;
rho(1) = 1.225;

%%velocity lift and drag coeff arrays for max excess power calculation.
Vel = linspace(10,600,1000);
lift_coeff = (2*W(1)./(rho(1).*Vel.^2.*S));
drag_coeff = Cd_0 + (lift_coeff.^2)./(pi().*e_prime.*AR);

%power required and available arrays based on velocity.
Preq = 0.5.*rho(1).*drag_coeff.*S.*Vel.^3;
Pav = Ta(1).*Vel;

%Max excess power and index for finding velocity.
[EP_max,index] = max(Pav - Preq);

%Initial vertical velocity.
V_v(1) = EP_max/W(1);

%Updating intial conditions for the for loop.
Pr(1) = Preq(index);
Pa(1) = Pav(index);
Tr(1) = Pr(1)/Vel(index);
Cl(1) = lift_coeff(index);
Cd(1) =  drag_coeff(index);
V(1) = Vel(index);
h(1) = 0;
fuel_used(1) = 0;
distance(1) = 0;
theta(1) = asin((Ta_sl - Tr(1))/Wmax);

%% Climb 

for i=2:length(time)
    
    %%density and thrust required calculations:
    [T,a,P,rho(i)]  = atmosisa(real(h(i-1)));
    Ta(i) = Ta_sl*(rho(i)/rho_sl)^0.9;
   
    %%Weight and fuel consumed:
    W(i) = W(i-1) - c_t*Ta(i-1)*dt; %updating weight from last
    fuel_used(i) = fuel_used(i-1) + c_t*Ta(i-1)*dt;

    
    %%Maximum Excess power and velocity calculations:
    V_array = linspace(30,1000,5000); %possible velocity array.
    
    cl = @(vel)(2.*W(i).*cos(theta(i-1))./(rho(i).*vel.^2.*S)); %%lift coefficent function.

    EP = @(vel)(Ta(i) - 0.5.*rho(i).*S.*(Cd_0 + (cl(vel).^2)./(pi().*e_prime.*AR)).*vel.^2 - W(i-1).*sin(theta(i-1))).*vel; %%excess power.
    
    [EP_max,index] = max(EP(V_array)); %%max excess power and index.
    
    V_v(i) = EP_max/W(i); %%vertical velocity calculation.
    
    V(i) = V_array(index); %% air speed.
    
    V_h(i) = V(i)*cos(theta(i-1)); %%horizontal velocity.
    
    %%updating altitude and distance:
    h(i) = h(i-1) + V_v(i)*dt;
    distance(i) = distance(i-1) + V_h(i)*dt;
    
    %%lift and drag coefficients
    Cl(i) = cl(V(i));
    Cd(i) = Cd_0 + (Cl(i)^2)/(pi()*e_prime*AR);
    
    %%thrust required, power required, and available:
    Tr(i) = 0.5*rho(i)*S*(Cd_0 + (Cl(i)^2)/(pi()*e_prime*AR))*V(i)^2 + W(i)*sin(theta(i-1));
    Pr(i) = Tr(i)*V(i);
    Pa(i) = Ta(i)*V(i);
    
    %%updating the climb angle:
    theta(i) = asin(V_v(i)/V(i));
    
    %%break condition
    if(h(i) >= 9000)
        index = find(h == h(i));
        fprintf("loop broken");
        
        break;
    end
    
end

%% Initial conditions for cruise portion: 
Ta_cruise = 0.9*Ta(index);

%%Cruise velocity calculation:
V_array = linspace(100,500,1000);

Cl_cruise = (2*Wmax)./(rho(index).*S.*V_array.^2); %%cruise CL array.

Cd_cruise = Cd_0 + (Cl_cruise.^2)./(pi().*e_prime.*AR); %%Cruise CD array.

Velocity = sqrt(2.*Ta_cruise./(rho(i).*Cd_cruise.*S)); %% setting drag = thrust to find max velocity:

%%cruise velocity calculation.
V_cruise = max(Velocity); 

%Updating the weight value:
W(index+1) = W(index) - c_t*Ta(index)*dt;

fprintf("%.2f\n",V_cruise);

%% Cruise:
for i = index + 1 :length(time)
    
    V(i) = 245; %%cruise velocity in m/s.
    
    %%updating density, weight and fuel used:
    rho(i) = rho(i-1); 
    W(i) = W(i-1) - c_t*Ta(i-1)*dt; 
    fuel_used(i) = fuel_used(i-1) + c_t*Ta(i-1)*dt;
    
    %%Lift and drag coeffs.
    Cl(i) = 2*W(i)/(rho(i)*S*V(i)^2);
    Cd(i) = Cd_0 + (Cl(i)^2)/(pi()*e_prime*AR);
    
    %%Thrust and power calculations:
    Tr(i) = 0.5*rho(i)*Cd(i)*S*V(i)^2;
    Ta(i) = Tr(i);
    
    Pr(i) = Tr(i)*V(i);
    Pa(i) = Ta(i)*V(i);
    
    %climb angle(should be zero), horizonal, and vertical velocities:
    theta = asin(Ta(i) - Ta(i))/W(i); %climb angle (should be zero);
    V_v(i) = V(i)*sin(theta);
    V_h(i) = V(i)*cos(theta);
    
    %%updating height and distance:
    h(i) = h(i-1) + V_v(i)*dt;
    distance(i) = distance(i-1) + V_h(i)*dt;
    
    %%break condition:
    if (distance(i) >= 2828e3)
        index = i;
        break;
    end
        
end

%% Intial decent

for i = index+1 :length(time)
    
    %%updating density, weight and fuel used:
    [T,a,P,rho(i)]  = atmosisa(real(h(i-1)));
    W(i) = W(i-1) - c_t*Tr(i-1)*dt;
    fuel_used(i) = fuel_used(i-1) + c_t*Ta(i-1)*dt;
    
    %Air velocity and decent angle:
    V(i) = 0.7*V_cruise;
    theta(i) = 3*pi()/180;
    
    %%lift and drag coeffs:
    Cl(i) = 0.15;
    Cd(i) = 0.05 + (Cl(i)^2)/(pi()*e_prime*AR);
    
    %%thrust required and available:
    Tr(i) = (0.5*rho(i)*Cd(i)*S*V(i)^2) - W(i)*sin(theta(i));
    Ta(i) = Tr(i);
    
    %%vertical and horizontal velocities:
    V_v(i) = V(i)*sin(theta(i));
    V_h(i) = V(i)*cos(theta(i));
    
    %%updating altitude and distance:
    h(i) = h(i-1) - V_v(i)*dt;
    distance(i) = distance(i-1) + V_h(i)*dt;
    
    %%break condition:
    if (h(i) <= 3000)
        dist = distance(i);
        index = find(distance == dist);
        break;
        
    end
    
end


%% Waiting Pattern.

%%initializing time count:
time_count = zeros(length(time),1);
time_count(index) = 0;

%%initializing turn angle:
turn_angle = zeros(length(time),1);
turn_angle(index) = 0;

for i=index+1:length(time)
    
    %%updating density:
    [T,a,P,rho(i)]  = atmosisa(real(h(i-1)));
    
    %%updating weight and fuel used:
    W(i) = W(i-1) - c_t*Tr(i)*dt;
    fuel_used(i) = fuel_used(i-1) + c_t*Ta(i-1)*dt;
    
    %Calculating CL for max CL/CD.
    V_wait = linspace(40,300);
    Cl_wait = 2.*W(i)./(rho(i).*S.*V_wait.^2);
    Cd_wait = Cd_0 + (Cl_wait.^2)./(pi().*e_prime.*AR());
    
    [Cl_CD_MAX,Cl_index] = max(Cl_wait./Cd_wait); %%maximizing CL/CD:
    
    %%Lift, drag coeffs and air velocity:
    Cl(i) = Cl_wait(Cl_index); 
    Cd(i) = Cd_0 + Cl(i)^2/(pi()*e_prime*AR);
    V(i) = V_wait(Cl_index);
    
    %thrust and power, available and required:
    Tr(i) = 0.5*rho(i)*Cd(i)*S*V(i)^2;
    Ta(i) = Tr(i);
    Pr(i) = Tr(i)*V(i);
    Pa(i) = Tr(i)*V(i);
    
    %vertical climb angle: (should be zero)
    theta(i) = asin((Ta(i) - 0.5*rho(i)*Cd(i)*S*V(i)^2)/W(i));
    
    %Turn rate and angle in turn.
    omega = 2*pi()/1200; 
    turn_angle(i) = turn_angle(i-1) + omega*dt;
    
    
    %updating vertical and horizontal velocities, altitude, and distance.
    V_h(i) = V(i)*cos(turn_angle(i)); %X component of velocity.
    
    V_v(i) = V(i)*sin(theta(i)); %vertical velocity (should be zero);
    
    %%Updating altitude and distance:
    h(i) = h(i-1) + V_v(i)*dt;
    distance(i) = distance(i-1) + V_h(i)*dt;
    
    %%updating time count:
    time_count(i) =  time_count(i-1) + dt;
    
    %%break condition:
    if time_count(i) >= 1200
        index = find(time_count == time_count(i));
        
        break;
        
    end
    
end

%% Final decent:

for i = index: length(time)
    
    %%updating density, weight and fuel used:
    [T,a,P,rho(i)]  = atmosisa(real(h(i-1)));
    W(i) = W(i-1) - c_t*Tr(i)*dt; %%updating weight.
    fuel_used(i) = fuel_used(i-1) + c_t*Ta(i-1)*dt;
    
    
    theta(i) = 3*pi()/180; %%decent angle
    
    %%Cd and Cl calculations.
    Cl(i) = 2.5;
    Cd(i) = 0.07 + (Cl(i)^2)/(pi()*e_prime*AR);
   
    %%velocity calculations:
    V_stall = sqrt(2*W(i)/(rho(i)*Cl(i)*S));
    V(i) = 1.2*V_stall;
    V_h(i) = V(i)*cos(theta(i));
    V_v(i) = -V(i)*sin(theta(i));
    
    %%Thrust and power, required and available.
    Tr(i) = 0.5*rho(i)*Cd(i)*S*V(i)^2 - W(i)*sin(theta(i));
    Ta(i) = Tr(i);
    
    Pr(i) = Tr(i)*V(i);
    Pa(i) = Ta(i)*V(i);
    
    %%updating height, distance 
    h(i) = h(i-1) + V_v(i)*dt;
    distance(i) = distance(i-1) + V_h(i)*dt;
    
    %%break condition
    if h(i)<= 0
        
        fprintf("%.4f/n",distance(i));
        
        break;
        
    end
    
end

%% Trimming arrays:

fuel_used_plot = fuel_used(1:i);
Cl_plot = Cl(1:i);
Cd_plot = Cd(1:i);

Pr_plot = Pr(1:i);
Pa_plot = Pa(1:i);

V_plot = V(1:i);

altitude = h(1:i);
distance_plot = distance(1:i);


%% Plotting the results:
figure;
plot(distance_plot,altitude);
xlabel("distance [m]");
ylabel("altitude [m]");
hold on;

hold off;

figure;
plot(distance_plot,V_plot);
xlabel("distance [m]");
ylabel("air velocity [m]");
hold off;

figure;
plot(distance_plot,Pr_plot,'b','DisplayName',"power required");
xlabel("distance [m]");
ylabel("Power required [m]");
hold on;
plot(distance_plot,Pa_plot,'r','DisplayName',"power available");
legend;
hold on;
hold off;

figure;
plot(distance_plot,Cd_plot);
xlabel("distance [m]");
ylabel("Drag coefficient");
hold on;
hold off;



figure;
plot(distance_plot,Cl_plot);
xlabel("distance [m]");
ylabel("Lift coefficient");
hold on;
hold off;

figure;
plot(distance_plot,fuel_used_plot);
xlabel("distance [m]");
ylabel("fuel consumed [N]");
hold on;
hold off;


%%



