Implementation of a parallel version of the Fpop algorithm. (Credit: sequential implementation <a href="https://github.com/vrunge/fpop">Vincent Runge</a>)

Add this two lines in your live R session before compiling the project:

```r
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp"))
```