clc,clear;
close all;
%% parameter set up
omega_e=7.292e-5;               %omega of the earth


R=6371000;                      %radius of the earth
Longitude(1)=(110+57/60)*pi/180;         %longitude of the launch site
Latitude(1)=(19+37/60)*pi/180;           %latitude of the launch site
h(1)=100;       %initial height

g=9.79;        %initial gravity

V(1,:) = [0,0,0];

t=0.01;         %update time of gyro output
k=20;           %update ratio of T and t
T=0.2;          %update time of rotation of the geographical frame in second
time=1460;         %data time in second
num = time/T;                      %7300

zero_angle = [90,0,-90];
qx=[cosd(zero_angle(1)/2) sind(zero_angle(1)/2) 0 0];
qy=[cosd(zero_angle(2)/2) 0 sind(zero_angle(2)/2) 0];
qz=[cosd(zero_angle(3)/2) 0 0 sind(zero_angle(3)/2)];
q(1,:)=quatmultiply(quatmultiply(qz,qx),qy);

load mission.mat;               %load the data
%% iteration
for count_1=1:num
    t(count_1)=count_1;
    for count_2=1:k
		w=0.01/(3600*57.296)*[GGM((count_1-1)*k+count_2,1),...
			GGM((count_1-1)*k+count_2,2),GGM((count_1-1)*k+count_2,3)];
		norm_w=norm(w);
		ome_w=[0 -w(1) -w(2) -w(3);
			  w(1) 0 w(3) -w(2);
			  w(2) -w(3) 0 w(1);
			  w(3) w(2) -w(1) 0];
		q(count_1,:)=(((1-norm_w^2/8+norm_w^4/384)*eye(4)+(0.5-norm_w^2/48)*ome_w)*q(count_1,:)')';
    end

    omega_rela=[-V(count_1,2)/(R+h(count_1)) V(count_1,1)/(R+h(count_1))+omega_e*cos(Latitude(count_1)) ...
        V(count_1,1)/(R+h(count_1))*tan(Latitude(count_1))+omega_e*sin(Latitude(count_1))];
    theta=T*omega_rela;
    tehta_re=norm(theta);
    ne=omega_rela/norm(omega_rela);
    q_g=[cos(tehta_re/2),sin(tehta_re/2)*ne];
 
    q(count_1+1,:)=quatmultiply(quatinv(q_g),q(count_1,:));   % quaternion updating
    angle_eular(count_1,:)=q2eul(q(count_1+1,:));        %transform and record angles
    q_f=quatmultiply(quatmultiply(q(count_1+1,:),[0 1e-7*g*AAM(count_1,:)]),...
        quatinv(q(count_1+1,:)));
    q_a=q_f(2:4)-cross(omega_rela+[0 omega_e*cos(Latitude(count_1)) ...
        omega_e*sin(Latitude(count_1))],V(count_1,:))-[0 0 g*R^2/(R+h(count_1))^2];  
    
    V(count_1+1,:)=V(count_1,:)+T*q_a;	%Updating velocity and location vectors 
    Longitude(count_1+1)=Longitude(count_1)+ ...
        T*0.5*(V(count_1,1)+V(count_1+1,1))/((R+h(count_1))*cos(Latitude(count_1)));
    Latitude(count_1+1)=Latitude(count_1)+T*0.5*(V(count_1,2)+V(count_1+1,2))/(R+h(count_1));
    h(count_1+1)=h(count_1)+T*0.5*(V(count_1,3)+V(count_1+1,3));
end
%% display data and show image
%Longtitude-latitude
figure(1)
title('Longtitude-latitude');
xlabel('Longitude/degree');
ylabel('latitude/degree');
hold on;
plot(Longitude/pi*180,Latitude/pi*180);
grid on;
%height-time 
figure(2)
xlabel('time/s')
ylabel('height/m')
title('height-time');
hold on;
plot([0:7300]/5,h);
grid on;
%Longitue and latitude of the rocket
display(Longitude(num+1),'Longitude');
display(Latitude(num+1),'Latitude');
display(h(num+1),'Height');
%display the final velocity
display('The velocity at the end of 1460s is');
display(V(num+1,:),'V_end');
%Quaternion after 1460s
Quaternion=q(num+1,:)
%show the track of rocket in a earth figure
open('globe.fig')
hold on
plot3m(Latitude/pi*180,Longitude/pi*180,h/R);
figure(4)
title('Euler angle');
xlabel('Time /s');
ylabel('Angle /degree');
grid on;
hold on;
plot([1:7300]/5,angle_eular()/pi*180);
