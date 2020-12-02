Implementation of a parallel version of the Fpop algorithm. 

Add these two lines in your live R session before compiling the project:

```r
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")
```
