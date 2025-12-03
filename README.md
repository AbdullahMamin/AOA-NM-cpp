## AOA-NM

This is an implementation of AOA-NM from the work of Betul, et al. Source: https://doi.org/10.1016/j.knosys.2023.110554

## How to run

Simply use the provided makefile to build the program. There are no dependencies, all that is needed is a C++ compiler that can compile the C++14 standard. You can then run the built executable "out" in the terminal and see the resulting CSV files after solving for the Beale function example. If you wish to optimize another function, look at main.cpp for instructions on how to do so.

## How to use

Refer to the comments in main.cpp to see how it is used. Essentially, you instantiate a solver object with some parameters and the objective function. Then you call solve to get your result. You can optionally get an evaluation metric if you wish to evaluate how many objective function calls it took to get to the result found.

## License

The code is licensed under GPLv3+.
