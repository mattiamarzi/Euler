In this project I solve the shock-tube problem in one dimension using high resolution schemes.

When we try to solve the system of Euler differential equations numerically, if we want a precision of the second order in the spatial discretization, we make an error which generates a dispersive effect around discontinuites.
We can see this effect when we try to use the Lax-Wendroff scheme:

![LW](https://github.com/mattiamarzi/Shock-Tube-Problem/assets/133958148/dfecb326-ad9a-4b7f-b231-8e37a3cc8d4d)

To solve this problem we can use a method with a spatial precision of the first order (such as the Roe scheme) and then increase the spatial precision to the second order (MUSCL scheme) using flux-limiters (such as the MINMOD and the SUPERBEE). In this way the dispersive effect is completely solved:

![M_MIN](https://github.com/mattiamarzi/Shock-Tube-Problem/assets/133958148/8523097e-6ce8-40be-a78e-b3149e99c6f6)
