# Known sources of errors
We summarize the major sources of errors so the user is aware of them:
- inserting points into the grid that are outside of the prescribed domain. This will cause in the best case weird verteces that do not exist or - in most cases - it will crash the algorithm. However, since the algorithm allows for unbounded domains, this is easy to prevent.
- The calculation of verteces and volumes is precise and reproduceable for regular grids. However, for nodes in non-general position, `VI_POLYGON` volume calculations tend to have a worst case error in the range of 2.0% in 5D and 6D and this value may increase in higher dimensions. You can test this on a unit cube to get an idea for a specific dimension. 
