# lowrankquadrmax
An algorithm to (approximately) solve sparse and nonnegative PCA
quadratic maximization problems

The algorithm operates on the raw data (samples x features matrix) or
the covariance matrix. 

It is based on Algorithm 3 in [Nonnegative Sparse PCA with Provable Guarantees](http://jmlr.csail.mit.edu/proceedings/papers/v32/asteris14.pdf) (ICML 2014).

See [QuickDemo.m](matlab/bin/QuickDemo.m) for an example on how to use the code.

### External code
The log4m.m created by Luke Winslow and is available through Matlab Central, [here](http://www.mathworks.com/matlabcentral/fileexchange/37701-log4m-a-powerful-and-simple-logger-for-matlab)).
The file is distributed under the BSD license.

