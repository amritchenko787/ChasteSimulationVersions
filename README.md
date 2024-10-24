1. To install dependencies for runnng Chaste, follow instructions here - https://chaste.github.io/docs/.
2. After downloading Chaste source code in the directory, verify installation by following instructions here - https://chaste.github.io/docs/user-guides/cmake-first-run/. 
3. To run any of the Chaste tutorials, it is assumed you have first configured Chaste using CMake. Follow instructions here - https://chaste.github.io/docs/user-tutorials/.
4. To run Durotaxis, Haptotaxis simulations follow steps - 
 i cd /path/to/chaste_build.
 ii make TestRunningVertexBasedSimulationTutorial_ab1.
 iii ctest -V -R TestRunningVertexBasedSimulationTutorial_ab1.
5. To visualize simulations in Paraview, follow instructions here - https://chaste.github.io/docs/user-tutorials/visualizingwithparaview/.

