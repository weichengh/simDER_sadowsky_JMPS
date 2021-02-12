To run to code, you may need a Linux system and install the following things:

1, eigen3, https://eigen.tuxfamily.org/dox-devel/index.html

2, openGL, sudo apt-get install freeglut3-dev

3, lapack, sudo apt-get install libatlas-base-dev

4, gfortran, sudo apt-get install gfortran

5, make:
g++ -I /usr/local/include/eigen-3.3.7/ main.cpp world.cpp elasticRod.cpp elasticStretchingForce.cpp elasticRibbonForce.cpp externalGravityForce.cpp inertialForce.cpp dampingForce.cpp timeStepper.cpp setInput.cpp -llapack -lGL -lglut -lGLU -Ofast -o simDER

6, run the code:
./simDER option.txt

To modify the boundary condition, you can go to "world.cpp" line 172-182, and line 188-234.

Feel free to email me if you have any questions.

