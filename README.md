# mf_secret_sharing
A mini-project to show that it is feasible to use modular forms instead of
polynomials for secret sharing.

## Usage

Make sure to install Sagemath (e. g. `#pacman -S sagemath` in Arch).

Then to run the test suite, open Sagemath in the directory of the repository and
run:

```
sage: attach("test_suite.sage")
sage: runtests_random_r()
```

You will most likely encounter the error
`ValueError; Invalid choice of r: [nbr]`. This is nothing to worry about, and just
indicates that the coefficient `r` in `create_public_data` satisfied $c(f_0;r)=0$.
All that needs to be done then is to change the `r`.

Why $c(f_0;r)=0$ is bad will be explained in an upcoming formal write-up.

## Todos

We can actually reconstruct the vote by using less than "Sturm bound" number of
coefficients. This means that Shamir's criterion

> but even complete knowledge of `k - 1` pieces reveals absolutely no
  information about `D` [ed: the vote sum]

is violated. However, we can just replace the Sturm bound with the dimension
instead. However, we need to prove that this is ok.
