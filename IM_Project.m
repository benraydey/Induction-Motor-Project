clear;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Energy Efficient (EE) Motor parameters')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
 
f=50;                    %Supply frequency [Hz]
p=4;                     %Number of poles
V1=380/sqrt(3);          %Supply voltage [phase]
R1_ee = 1.5;             %Stator winding resistance [ohms/phase]
X1_ee = 3.642;           %Stator winding leakage reactance [ohms/phase]
Xm_ee = 72.252;          %Stator winding magnetising reactance [ohms/phase]
X2p_ee = 3.642;          %Rotor winding leakage reactance reffered to stator [ohms/phase]
R2p_ee = 1.994;          %Rotor winding resistance reffered to stator [ohms/phase]

fprintf('f=%f\n',f);
fprintf('p=%f\n',p);
fprintf('V1=%f\n',V1);
fprintf('R1_ee=%f\n',R1_ee);
fprintf('X1_ee=%f\n',X1_ee);
fprintf('Xm_ee=%f\n',Xm_ee);
fprintf('X2p_ee=%f\n',X2p_ee);
fprintf('R2p_ee=%f\n',R2p_ee);

%%
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Standard Efficiency (SE) Motor parameters')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

R1_se = 2.087;           %Stator winding resistance [ohms/phase]
X1_se = 4.274;           %Stator winding leakage reactance [ohms/phase]
Xm_se = 66.560;          %Stator winding magnetising reactance [ohms/phase]
X2p_se = 4.274;          %Rotor winding leakage reactance reffered to stator [ohms/phase]
R2p_se = 2.122;          %Rotor winding resistance reffered to stator [ohms/phase]

fprintf('f=%f\n',f);
fprintf('p=%f\n',p);
fprintf('V1=%f\n',V1);
fprintf('R1_se=%f\n',R1_se);
fprintf('X1_se=%f\n',X1_se);
fprintf('Xm_se=%f\n',Xm_se);
fprintf('X2p_se=%f\n',X2p_se);
fprintf('R2p_se=%f\n',R2p_se);
%%
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 1: Thevenin Equiv Cct Parameters for EE Motor:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

Vth_ee=Xm_ee/sqrt(R1_ee^2+(X1_ee+Xm_ee)^2)*V1;         %Thevenin equiv voltage source [V] (Equ 5.45 - Sen)
Zth_ee=1i*Xm_ee*(R1_ee+1i*X1_ee)/(R1_ee+1i*(X1_ee+Xm_ee));   %Thevenin equiv impedance
Rth_ee=real(Zth_ee);                          %Thevenin equiv resistance [ohms]
Xth_ee=imag(Zth_ee);                          %Thevenin equiv reactance [ohms]


fprintf('Vth=%f\n',Vth_ee);
fprintf('Rth=%f\n',Rth_ee);
fprintf('Xth=%f\n',Xth_ee);

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Question 1: Thevenin Equiv Cct Parameters for SE Motor:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

Vth_se=Xm_se/sqrt(R1_se^2+(X1_se+Xm_se)^2)*V1;         %Thevenin equiv voltage source [V] (Equ 5.45 - Sen)
Zth_se=1i*Xm_se*(R1_se+1i*X1_ee)/(R1_se+1i*(X1_se+Xm_se));   %Thevenin equiv impedance
Rth_se=real(Zth_se);                          %Thevenin equiv resistance [ohms]
Xth_se=imag(Zth_se);                          %Thevenin equiv reactance [ohms]


fprintf('Vth=%f\n',Vth_se);
fprintf('Rth=%f\n',Rth_se);
fprintf('Xth=%f\n',Xth_se);

%%
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('QUESTION 2:Torque versus speed characteristics for EE and SE Motors:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

ns=120*f/p;         %Synchronous speed [rpm]
ws=2*pi*ns/60;      %Synchronous speed [rad/sec]
s=0.0005:0.0005:1;  %Slip [pu]
n=(1-s)*ns;         %Rotor speed [rpm]
w=2*pi*n/60;        %Rotor speed [rad/sec]
 

% Torque calculation for EE motor
Tmech_ee = 3/ws*Vth_ee^2./((Rth_ee+R2p_ee./s).^2+(Xth_ee+X2p_ee)^2).*R2p_ee./s;    %Total Tmech = {3*(Equ5.54 - Sen)}

% Torque calculation for SE motor
Tmech_se = 3/ws*Vth_se^2./((Rth_se+R2p_se./s).^2+(Xth_se+X2p_se)^2).*R2p_se./s;

%Plot the Torque vs speed characteristics of both motors on same plot:
%subplot(2,2,1),
plot(n,Tmech_ee),xlabel('n [rpm]'),ylabel('Torque [Nm]'),...
    title('Torque vs speed'),grid on,...
    hold on
plot(n,Tmech_se)

hold off

%%
% Calculate starting torque of machines
disp('(a) Starting Torque')

% Starting torque of EE motor
Tstart_ee = 3/ws*Vth_ee^2/((Rth_ee+R2p_ee/1)^2+(Xth_ee+X2p_ee)^2)*R2p_ee/1;    % Tstart = Tmech @ s = 1

% Starting torque of SE motor
Tstart_se = 3/ws*Vth_se^2/((Rth_se+R2p_se/1)^2+(Xth_se+X2p_se)^2)*R2p_se/1;

fprintf('Starting torque of EE motor = %f N.m\n',Tstart_ee);
fprintf('Starting torque of SE motor = %f N.m\n',Tstart_se);

disp(['The starting torque can be increased by increasing the rotor resistance, and' ...
    ' decreased by decreasing the rotor resistance.'])
%%
% Calculate maximum (breakdown) torques of machines
disp('(b) Maximum Torque')

% Maximum torque of EE motor
Tmax_ee = 3/(2*ws)*Vth_ee^2/( Rth_ee + sqrt(Rth_ee^2 + (Xth_ee+X2p_ee)^2) )

% Maximum torque of SE motor
Tmax_se = 3/(2*ws)*Vth_se^2/( Rth_se + sqrt(Rth_se^2 + (Xth_se+X2p_se)^2) )

fprintf('Maximum torque of EE motor = %f N.m\n',Tstart_ee);
fprintf('Maximum torque of SE motor = %f N.m\n',Tstart_se);

disp(['The maximum torque is dependent of the the Thevenin equivalent voltage Vth. Thus, the maximum torque can ' ...
    'be increased by increasing Vth and decreased by decreasing Vth'])
%%
% Caculate speed (in rpm) at which maximum torque occurs for the machines
disp('(c) Speed at which maximum torque occurs')

% EE motor
s_Tmax_ee = R2p_ee/sqrt(Rth_ee^2 + (Xth_ee + X2p_ee)^2)
n_Tmax_ee = ns*(1-s_Tmax_ee)

% SE motor
s_Tmax_se = R2p_se/sqrt(Rth_se^2 + (Xth_se + X2p_se)^2)
n_Tmax_se = ns*(1-s_Tmax_se)

fprintf('Speed at maximum torque of EE motor = %f rpm\n',n_Tmax_ee);
fprintf('Speed at maximum torque of SE motor = %f rpm\n',n_Tmax_se);


%%
% Stator current vs. speed characteristics of machines
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('QUESTION 3:Stator current versus speed characteristicc of EE and SE Motors:')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


Tmech_se = 3/ws*Vth_se^2./((Rth_se+R2p_se./s).^2+(Xth_se+X2p_se)^2).*R2p_se./s;

% EE motor
Z1_ee = R1_ee + 1i*X1_ee + 1i*Xm_ee*(R2p_ee./s+1i*X2p_ee)./(R2p_ee./s+1i*(Xm_ee+X2p_ee))
I1_ee = V1./Z1_ee

% SE motor
Z1_se = R1_se +  1i*X1_se + 1i*Xm_se*(R2p_se./s+1i*X2p_se)./(R2p_se./s+1i*(Xm_se+X2p_se))
I1_se = V1./Z1_se

% Plot stator current vs. speed characteristic curves
plot(n,abs(I1_ee)),xlabel('n [rpm]'),ylabel('Stator current [A]'),...
    title('Stator current vs speed'),grid on,...
    hold on
plot(n,abs(I1_se))

hold off
%%
% (a) Stator current at startup
disp('(a) Stator current at startup')

% EE motor
fprintf('Stator current of Energy Efficient machine at startup = %f A',abs(I1_ee(2000)))

% SE motor
fprintf('Stator current of Standard machine at startup = %f A',abs(I1_se(2000)))
%%
% b) Stator current - no load and full load conditions
disp('(a) Explain how the stator current at start-up would change when the machine is started under no-load and full-load conditions?')
disp(['Stator current would be at a maximum under full load condtions and at a minimum and no load conditions.'])
%%
% (c) Stator currents when the machines develop maximum torque
disp('(c) Stator currents when the machines develop maximum torque')

% EE motor
Z1_sTmax_ee = R1_ee + 1i*X1_ee + 1i*Xm_ee*(R2p_ee./s_Tmax_ee+1i*X2p_ee)./(R2p_ee./s_Tmax_ee+1i*(Xm_ee+X2p_ee))
I1_sTmax_ee = V1./Z1_sTmax_ee
fprintf('Stator current of Energy Efficient machine at max torque = %f A',abs(I1_sTmax_ee))

% SE motor
Z1_sTmax_se = R1_se + 1i*X1_se + 1i*Xm_se*(R2p_se./s_Tmax_se+1i*X2p_se)./(R2p_se./s_Tmax_se+1i*(Xm_se+X2p_se))
I1_sTmax_se = V1./Z1_sTmax_se
fprintf('Stator current of Standard machine at max torque = %f A',abs(I1_sTmax_se))
%%
% (d) Stator current when the machines are operating under no-load conditions
disp('(d) Stator current when the machines are operating under no-load conditions')


%%
% 4. Power factor vs. speed characteristic of machines

% EE motor
PF_ee = cos(angle(I1_ee))

% SE motor
PF_se = cos(angle(I1_se))

% Plot power factor vs. speed characteristic curves
plot(n,PF_ee),xlabel('n [rpm]'),ylabel('Power factor'),...
    title('Power factor vs speed'),grid on,...
    hold on
plot(n,abs(PF_se))

hold off
%%
% (a) Power factors at startup

PF_startup_ee = PF_ee(2000)
PF_startup_se = PF_se(2000)

fprintf('Power factor at startup of EE motor = %f \n',PF_startup_ee);
fprintf('Power factor at startup of SE motor = %f \n',PF_startup_se);

%%
% (b) Power factors when the machines develop maximum torque



PF_Tmax_ee = cos(angle(I1_sTmax_ee))
PF_Tmax_se = cos(angle(I1_sTmax_se))

fprintf('Power factor at maximum torque of EE motor = %f \n',PF_Tmax_ee);
fprintf('Power factor at maximum torque of SE motor = %f \n',PF_Tmax_se);
%%
% (c) Power factor when the machines are operating under no-load conditions

PF_noload_ee = PF_ee(1)
PF_noload_se = PF_se(1)

fprintf('Power factor at no load of EE motor = %f \n',PF_noload_ee);
fprintf('Power factor at no load of SE motor = %f \n',PF_noload_se);

%%
% (d) Determine the best power factors that these machines can operate at

PF_best_ee = max(PF_ee)
PF_best_se = max(PF_se)

fprintf('Best power factor for EE motor = %f \n',PF_best_ee);
fprintf('Best power factor for SE motor = %f \n',PF_best_se);
%%
% (e) Determine the speeds (in rpm) at which the best power factors occur

pos_PF_best_ee = find(PF_ee == max(PF_ee))
n_best_ee = n(pos_PF_best_ee)

pos_PF_best_se = find(PF_se == max(PF_se))
n_best_se = n(pos_PF_best_se)

fprintf('Speed at which best power factor occurs for EE motor = %f \n',n_best_ee);
fprintf('Speed at which best power factor occurs for SE motor = %f \n',n_best_se);

%%



%% 
%