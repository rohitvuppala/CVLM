using DifferentialEquations
using Plots

# State Vector
#  X = [u,v,w,p,q,r,phi,theta,psi,x,y,z]

# time span
t_begin = 0
t_end = 10.0

tspan = (t_begin,t_end)

##########################################
# initial conditions
##########################################
u0 = 100
v0 = 0
w0 = 3 

p0 = 0
q0 = 0
r0 = 0 

phi0 = 0
theta0 = deg2rad(2)
psi0 = 0

x0 = 0
y0 = 0
z0 = -500


###############
# Let's get it running with CVLM in the loop - then, we will need to come up with a trim routine here for steady level flight. I've done this before, so I can add it. 
##############

# inital condition vector [tranlational velocities | angular velocities | euler angles | cartesian position ]
y0 = [u0;v0;w0;p0;q0;r0;phi0;theta0;psi0;x0;y0;z0]

#############################################
# parameter definitions - 
## **** Rohit - Can you make these into a Tuple? I tried unsuccessfully *******
#############################################
m = 10
g = 9.81
rho = 1.225
T = 0.2 # thrust (will need to determine from induced drag at trim?)
# Inertia constants for calculating GAMMA, C1-C9
Ix = 1
Iy = 1
Iz = 1
Ixz = 0.2
# constants to deal with inertial coupling
GAMMA = Ix*Iz-Ixz^2
c1 = [(Iy-Iz)*Iz-Ixz^2]/GAMMA
c2 = [(Ix-Iy+Iz)*Ixz]/GAMMA
c3 =  Iz/GAMMA
c4 = Ixz/GAMMA
c5 = (Iz-Ix)/Iy
c6 = Ixz/Iy
c7 = 1/Iy
c8 = [(Ix-Iy)*Ix-Ixz^2]/GAMMA
c9 = Ix/GAMMA
# *** Rohit -  add matrices from CVLM here
# all we need to do is define Px, Py, Pz, Qx, Qy, Qz as parameter to ODE solver
Px = zeros(6,6)
Py = zeros(6,6)
Pz = zeros(6,6)
Qx = zeros(6,6)
Qy = zeros(6,6)
Qz = zeros(6,6)

# rotating mass along longitudinal axis (spinning propellor or turbomachinery) - these will probably always be zero for using
Ip = 0  # rotating inertia 
OmegaP = 0  # rpm of rotating machinery



##################################################
# Define parameter vector (would be more convienient as a tuple)
##################################################
p = (m,g,rho,c1,c2,c3,c4,c5,c6,c7,c8,c9,T,Ip,OmegaP)


## ODE Function - For some reason I can't save this as a separate file.... probably Julia syntax 
function calc_dx(dx,x,p,t)
    
    #X = [u,v,w,p,q,r,phi,theta,psi,x,y,z]
    
    # parameters: 
    m,g,rho,c1,c2,c3,c4,c5,c6,c7,c8,T,Ip,OmegaP = p
    # current state values: 
    u,v,w,p,q,r,phi,theta,psi,x,y,z = x

    # Aerodynamic Force/Moment Calculation (without controls {for now})
    X = [u;v;w;p;q;r]

    Fx = rho*X'*Px*X
    Fy = rho*X'*Py*X
    Fz = rho*X'*Pz*X

    Mx = rho*X'*Qx*X
    Mx = rho*X'*Qy*X
    Mz = rho*X'*Qz*X

    
    # state derivative calculations

    dx[1] = r*v-q*w+Fx/m-g*sin(theta)+T/m
    dx[2] = p*w-r*u+Fy/m+g*cos(theta)*sin(phi)  # Force Equations
    dx[3] = q*u-p*v+Fz/m+g*cos(theta)*cos(phi)

    dx[4] = (c1*r+c2*p-c4*Ip*OmegaP)*q.+c3*Mx.+c4*Mz
    dx[5] = (c5*p+c7*Ip*OmegaP)*r-c6*(p^2-r^2)+c7*My   # moment equations
    dx[6] = (c8*p-c2*r-c9*Ip*OmegaP)*q+c9*Mz+c4*Mx

    dx[7] = p+tan(theta)*(q*sin(phi)+r*cos(phi))
    dx[8] = q*cos(phi)-r*sin(phi)
    dx[9] = (q*sin(phi)+r*cos(phi))/cos(theta)

    dx[10:12] = [1 0 0; 0 cos(phi) sin(phi);0 -sin(phi) cos(phi)]*[cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)]*[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1]*x[1:3]


end


prob = ODEProblem(calc_dx,y0,tspan,p)
sol = solve(prob)

plot(sol)





