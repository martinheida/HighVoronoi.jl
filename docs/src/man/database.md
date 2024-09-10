# Voronoi: Database Structure

There are currently three types of internal data storage implemented. They may be called using the keyword `vertex_storage` at the `VoronoiGeometry` constructor.

## Standard solution

The most recent and most efficient method is `DatabaseVertexStorage()`. It stores all information in one centralized database and uses a sophisticated indexing system for fast access to any information from all points in the code. It might be slightly slower than the other two methods for large grids in high dimension, but it requires much less memory and for smaller grids it is even faster than the other two methods. 

This is the only database that is reliably compatible with multithreading and as such it will be automatically enforced once `threading=MutliThread(...)` is provided as an option!

## Deprecated solution

Another option is the `ReferencedVertexStorage()` which is slower but may be useful in low dimensions. It has a decentralized dictionary-based data structure with additional references. It was first implemented to save memory compared to the very initial solution. 

## Initial solution

The `ClassicVertexStorage()` which is fast for integration algorithms in low dimensions and which was the first database structure underlying the computations. It builds solely on a system of dictionaries. For some applications it is efficient but it requires a lot of memory which becomes troublesome in high dimensions. However, for documentation of the development of the library, it is still useable.

