# Auto-tune-Wavenet-PID-Tail-Sitter-Control
Controlling Tail Sitter UAV by using radial base neural networks and IIR digital filters in order to perform system identification to tunne PID controllers

- The file named Tail_sitter_aut_2016a.slx works as a main;
- WPSO_PID2 is a function were identification and control algorithm were coding;
- motor.m contains thrust and torque polinomial fit of each motor;
- Aero_coefs2 contains aerodynamics coefficients of Tail-sitter aircraft for several angles of attack and Sideslip
- order_decent makes sort over PSO algorithm population after each time sample
