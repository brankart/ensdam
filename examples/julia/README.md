## Example of Julia interface to EnsDAM

### Warning: This interface is still experimental.

The most appropriate interface is based on the C-callable wrappers of the Fortran modules,
as illustrated in the examples:

- example_scores.jl
- example_probability.jl

The code can also be used in Julia by applying driectly the Python interface,
but this does not work for all functions. See example in:

- example_random.jl
