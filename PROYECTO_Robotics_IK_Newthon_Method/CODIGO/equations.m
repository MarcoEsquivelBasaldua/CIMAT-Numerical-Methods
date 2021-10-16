syms l1 l2 l3 theta1 theta2 theta3 a1

A1 = compute_dh_matrix(0,pi/2,l1,theta1)
A2 = compute_dh_matrix(l2,0,0,theta2)
A3 = compute_dh_matrix(l3,0,0,theta3)

A1*A2*A3