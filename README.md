# Spin-Up-Geomtry
This is a MATLAB software that solves the ferrohydrodynamic problem, using finite difference equations. 
The software does not require simplifications that limit its validity. So, it can be used for different field intensities and frequencies.

In the Matlab prompt, this must be given to the software in the following order: 

  SpinUpFlowDescriptor(nder,ndeθ,ndt ,tf,maxIterT,varargin)

where varargin relates to two optional inputs, considered as follows: 

a. If no argument is given (i.e. length(varargin) = 0), the software carries out a    
   simulation with default parameters. Such a scenario considers magnetic field    
   densities BmT = [3.4 4.5 5.6 6.8 7.9 9.0 10.1 11.3] mT, and the properties of      
   ferrofluid WBF1 (see Table 1 of the article),

b. In the case of a single input, it must be a vector of magnetic field densities     
   (in mT) for running the simulation,

c. Should there be two arguments, the first one must relate to the aforementioned    
   vector, whilst the other one must be a structure with the parameters of the    
   ferrofluid that will be simulated, and

d. Finally, in the case there are three or more arguments, only the first two are    
   considered and must comply with the requirements described in the previous      
   item.

