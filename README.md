
```
                    _____                ________________                   
                   |_   _|        /\    / ______________/                   
                     | |  __ _   /  \  | (___ |  ____/                      
                     | | / _' | / /\ \  \___ \| |__                         
                    _| || (_| |/ ____ \     ) |  __|                        
                   |_____\__, /_/   _\_\____) | |                           
                          __/ |   /__________/|_|                           
                         |___/                                              
                                                                            
 -------------------------------------------------------------------------- 
                Isogeometric-Assembling-by-Sum-Factorization                
 -------------------------------------------------------------------------- 
                                                                            
 https://github.com/IgASF/IgASF                                             
                                                                            
 Type                                                                       
                                                                            
   make all               ... to make executables of this library (debug)   
   make all-release       ... to make executables of this library (release) 
   make so                ... to make this library as so-file (debug)       
   make so-release        ... to make this library as so-file (release)     
   make unittest          ... to make and run unittests (debug)             
   make unittest-release  ... to make and run unittests (release)           
   make clean             ... to clean the build environment (debug)        
   make clean-release     ... to clean the build environment (release)      
   make echoall           ... to get the contents of the given variables    
                                                                            
 Use the files                                                              
                                                                            
   config/debug   and   config/release                                      
                                                                            
 to specify the corresponding configuration. If those files are deleted,    
 they will be restored based on the corresponding templates.                
                                                                            
 Invoke after 'make all'                                                    
                                                                            
   ./generateTest-debug   ... to create a configuration file for assembling 
                              the problem of interst                        
                                                                            
   ./runTest-debug test   ... to assemble a problem specified by such a     
                              configuration file test, like with            
                                                                            
   ./runTest-debug tests/QuarterAnnulus_stiff_128/deg03 -o mat              
                                                                            
                              to assemble and store result in file mat using
                              a binary format                               
                                                                            
   ./echoMatrix-debug mat ... to echo the matrix data from file mat         
                                                                            
   ./compareMatrices-debug mat1 mat2                                        
                          ... to compare such matrix files                  
                                                                            
   ./runTestsAndCompare-debug test[s]                                       
                          ... to assemble the problem with all methods and  
                              compare the results and write results to      
                              log.txt and compare.txt, like with            
                                                                            
   ./runTestsAndCompare-debug tests/QuarterAnnulus_stiff_128/*              

```

The code is based on the paper

A. Bressan, S. Takacs. *Sum-factorization techniques in Isogeometric Analysis*.
Submitted. [arXiv:1809.05471](https://arxiv.org/abs/1809.05471)

The first author has received funding from the European Research Council under
the European Union's Seventh Framework Programme (FP7/2007-2013) / ERC grant
agreement 339643. The second author was partially funded by the Austrian Science
Fund (FWF) under grant NFN S117-03.


