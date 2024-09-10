# [Multithreading](@id multithreading)

you can use the following keywords either in `threading` in the `search_settings` to occasionally use multi threading in Voronoi computations or you can set it for `integrate=....`. Note that `integrate=true` is equivalent to `integrate=SingleThread()`

## Periodic Meshes

Mutlithreading is currently not available for `fast=true` generation of periodic meshes!! 

## Automatic Inference of Threads

`threading=AutoThread()` in the `search_settings` will automatically infer the maximal number of availabler threads and correspondingly use `MultiThread` or `SingleThread` from below.


## Single Threaded Computations (Even if Julia is started with more threads)

Single threaded calculations can be enforced with `threading=SingleThread()` in the `search_settings`. Even if you start Julia with several threads but want to do a single threaded computation, it is strongly advised to use this option as this will call a specialized version of the code.

## Multi Threaded Computations

Parallelized compuations can be enforced with `threading=MultiThread(a,b)` in the `search_settings`. `a` and `b` provide parallelization information on an "outer" and an "inner" parallelization. While the use of `a` is save, it is currently not adviced to use `b>1`. If `a>Threads.nthreads()` it will be internally reduced to the maximally available number of threads. 

Currently, `b` is implemented as a parameter, because it is technically doable. However, more tests are needed before it can be said if it is advantageous compared to the sole usage of `a`.
