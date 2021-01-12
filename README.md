# Dynamic-scaling

This is the MATLAB code for the paper to appear in the Journal of Experimental Biology: 
Estimation of sinking velocities using free-falling dynamically scaled models: foraminifera as a test case
by Matthew Walker, Jörg U. Hammel, Fabian Wilde, Tatjana Hoehfurtner, Stuart Humphries and Rudi Schuech

The main script is Dynamic_Scaling.m which calls the other functions within this folder. The primary purpose of the code is to predict a scale at which to 3D print an enlarged model to replicate the dynamics of a real organism or particle sinking due to gravity. This is done either using the known Cd(Re) relationship for a sphere or by using a spline fit to existing experimental data collected thus far for enlarged model particles. Existing data needs to include measured sinking velocities, as well as basic parameters of the particle (i.e. length parallel to the flow, projected area and volume of material) and the fluid (density, viscosity). All values should be given in SI units.  The code also predicts the sinking velocity of the real particle, given the experimentally observed sinking velocities of a number of scaled models.

In the first instance contact either Matthew Walker (matthewwalkerbio@gmail.com) or Rudi Schuech (rudi.schuech@gmail.com) for assistance. Examplar data avaliable on request.
