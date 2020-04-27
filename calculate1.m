function R0_test = calculate1(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)
R0_test = alfa/r1+beta*epsilon/(r1*r2)+gamma*zeta/(r1*r3)+delta*eta*epsilon/(r1*r2*r4)+delta*zeta*theta/(r1*r3*r4)
end