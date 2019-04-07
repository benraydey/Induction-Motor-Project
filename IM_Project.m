clear;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Energy Efficient (EE) Motor parameters')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
 
f=60;                    %Supply frequency [Hz]
p=4;                     %Number of poles
V1=240/sqrt(3);          %Supply voltage [phase]
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
disp('QUESTION 2:Torque versus speed characteristics for EE Motor:')
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
disp('(c) Speed at maximum torque');

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
% (c) Stator currents when the machines develop maximum torque
disp('(c) Stator currents when the machines develop maximum torque')

% EE motor
Z1_sTmax_ee = R1_ee + 1i*X1_ee + 1i*Xm_ee*(R2p_ee/s_Tmax_ee+1i*X2p_ee)/(R2p_ee/s_Tmax_ee+1i*(Xm_ee+X2p_ee))
I1_sTmax_ee = V1/Z1_sTmax_ee
fprintf('Stator current of Energy Efficient machine at max torque = %f A',abs(I1_sTmax_ee))

% SE motor
Z1_sTmax_se = R1_se + 1i*X1_se + 1i*Xm_se*(R2p_se/s_Tmax_se+1i*X2p_se)/(R2p_se/s_Tmax_se+1i*(Xm_se+X2p_se))
I1_sTmax_se = V1/Z1_sTmax_se
fprintf('Stator current of Standard machine at max torque = %f A',abs(I1_sTmax_se))
%%
% (d) Stator current when the machines are operating under no-load conditions

% when machine has no load, it will operate at synchronous speed. Thus, the
% no-load stator current can be calculated as the stator current when slip
% s = 0

disp('(d) Stator current when the machines are operating under no-load conditions')

% EE motor
I1_NL_ee = abs(I1_ee(1));
fprintf('Stator current of Energy Efficient machine operating under no-load conditions = %f A',I1_NL_ee);

% SE motor
I1_NL_se = abs(I1_se(1));
fprintf('Stator current of Standard machine operating under no-load conditions = %f A',I1_NL_se);
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
%%
% (b) Power factors when the machines develop maximum torque
%%
% (c) Power factor when the machines are operating under no-load conditions
%%
% (d) Determine the best power factors that these machines can operate at
%%
% (e) Determine the speeds (in rpm) at which the best power factors occur
%%
% 5. Plot total input power,
%               stator copper losses,
%               air gap power,
%               rotor copper losses and
%               shaft power
%                               vs speed, on the same set of axis
% *Rotational losses of machines are negligable
% *Seperate plots for each machine

disp('5. Power and losses plots of machines');

% EE motor
I2p_ee = Vth_ee./(Rth_ee+R2p_ee./s+1i*(Xth_ee+X2p_ee));      % Rotor current referred to primary

P1cu_ee = 3*I1_ee.^2*R1_ee;         % stator copper loss
Pag_ee = 3*I2p_ee.^2*R2p_ee./s;     % air gap power
P2cu_ee = 3*I2p_ee.^2*R2p_ee;       % rotor copper loss
Pshaft_ee = (1-s).*Pag_ee;          % shaft power Pshaft = Pmech since Prot is negligable
%Pin_ee = Pshaft_ee+P1cu_ee+P2cu_ee;            % total input power = Pag + P1cu
Pin_ee = 3*V1*I1_ee;

% generate plots of powers and losses
plot(n,abs(Pin_ee),n,abs(Pag_ee),n,abs(Pshaft_ee),n,abs(P1cu_ee),n,abs(P2cu_ee),'linewidth',2),xlabel('n [rpm]'),ylabel('Power [W]'),...
    title('Total input power, stator copper loss, air gap power, rotor copper loss and shaft power vs speed for Energy Efficient motor'),grid on,...
legend({'Total input power','Air gap power','Shaft power','Stator copper loss','Rotor copper loss'},'Location','best');


% SE motor
I2p_se = Vth_se./(Rth_se+R2p_se./s+1i*(Xth_se+X2p_se));      % Rotor current referred to primary

P1cu_se = 3*I1_se.^2*R1_se;         % stator copper loss
Pag_se = 3*I2p_se.^2*R2p_se./s;     % air gap power
P2cu_se = 3*I2p_se.^2*R2p_se;       % rotor copper loss
Pshaft_se = (1-s).*Pag_se;          % shaft power Pshaft = Pmech since Prot is negligable
%Pin_se = Pag_se+P1cu_se;            % total input power = Pag + P1cu
Pin_se = 3*V1*I1_se;

% generate plots of powers and losses
plot(n,abs(Pin_se),n,abs(Pag_se),n,abs(Pshaft_se),n,abs(P1cu_se),n,abs(P2cu_se),'linewidth',2),xlabel('n [rpm]'),ylabel('Power [W]'),...
    title('Total input power, stator copper loss, air gap power, rotor copper loss and shaft power vs speed for Standard motor'),grid on,...
legend({'Total input power','Air gap power','Shaft power','Stator copper loss','Rotor copper loss'},'Location','best');
%%
% 5. (a) Stator and rotor copper losses at start up
disp('(a) Stator and rotor copper losses at start up');

% EE motor
fprintf('Stator copper loss of Energy Efficient motor at startup = %f W',abs(P1cu_ee(2000)));
fprintf('Rotor copper loss of Energy Efficient motor at startup = %f W',abs(P2cu_ee(2000)));

% SE motor
fprintf('Stator copper loss of Standard motor at startup = %f W',abs(P1cu_se(2000)));
fprintf('Rotor copper loss of Standard motor at startup = %f W',abs(P2cu_se(2000)));
%%
% 5. (b) Stator and rotor copper losses under no-load conditions
disp('(b) Stator and rotor copper losses under no-load conditions');

% EE motor
fprintf('Stator copper loss of Energy Efficient motor under no-load conditions = %f W',abs(P1cu_ee(1)));
fprintf('Rotor copper loss of Energy Efficient motor under no-load conditions = %f W',abs(P2cu_ee(1)));

% SE motor
fprintf('Stator copper loss of Standard motor under no-load conditions = %f W',abs(P1cu_se(1)));
fprintf('Rotor copper loss of Standard motor under no-load conditions = %f W',abs(P2cu_se(1)));
%%
% 6. Efficiency vs. speed characteristics
disp('6. Efficiency vs. speed characteristics')

% EE motor
eff_ee = Pshaft_ee./Pin_ee;      % efficiency = Pout/Pin

% SE motor
eff_se = Pshaft_se./Pin_se;      % efficiency = Pout/Pin

% Plot the characteristics
plot(n,abs(eff_ee),n,abs(eff_se),'linewidth',2),xlabel('n [rpm]'),ylabel('Efficiency'),...
    title('Efficiency vs speed'),grid on,...
    %legend({'EE motor','SE motor'},'Location','best');
%%
% 6. (a) Efficiencies of machines at maximum torque
disp('(a) Efficiencies of the motors at maximum torque')

%%
% 6. (b) Maximum efficiency of machines
disp('(b) Maximum efficiency of machine')

%fprintf('Maximum efficiency of Energy Efficient motor = %f',)
%% 
%