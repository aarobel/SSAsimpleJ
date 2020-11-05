# SSAsimpleJ
Model for 1D marine ice sheet evolution (SSA approximation) in Julia, using numerical approach from Schoof 2007

Finite difference is the same as MATLAB impementation (https://github.com/aarobel/SSAsimpleM), except that Jacobian is derived through automatic differentiation instead of finite difference, which makes the solver more accurate and 4x faster than the equivalent MATLAB version.

Other implementation details coming soon...
