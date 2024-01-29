To compile:
g++ -o main 4401-81.cpp  otrMathModel.cpp main.cpp 

In main.cpp fill missile parametrs, initial motion parameters, stabilization system parametrs, target position(if need)

To get the affected area:
Fill all parameters without target position, compile and run.

To get parametrs of trajectory at specific target:
Fill all parameters. Uncomment ``cout`` in function ``update()``, ``get_res()`` in ``otrMathModel.cpp``. Uncomment ``la.get_res(dt, iter)`` and comment all after this in ``main.cpp``. Compile and run.   
