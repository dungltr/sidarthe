function R0_test = calculate2(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)
R0_test = (alfa*r2*r3*r4+epsilon*beta*r3*r4+gamma*zeta*r2*r4+delta*eta*epsilon*r3+delta*zeta*theta*r2)/(r1*r2*r3*r4)
end