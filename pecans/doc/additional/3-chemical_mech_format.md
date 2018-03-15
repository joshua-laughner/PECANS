# Format for mechanism input files

## The basics

Specifying a chemical mechanism requires at least two files:

- a `.spc` file that lists the species used in the reaction. (This helps catch typos in the next file.)
- a `.eqn` file that lists the chemical reactions that occur and their rate constants.

By requiring the species file, the mechanism generator can verify that all the species used in the
reactions are expected. If there is a typo, say "ISPO" instead of "ISOP" for isoprene, the parser
can alert you to this, because "ISPO" is not a defined species.

There is one additional file that can be provided for each mechanism, which is a `.rate` file, which 
provides custom kinetic rate constants for that mechanism (more information below).

When you execute `build_pecans.py`, it finds the corresponding `.spc`, `.eqn`, and `.rate` (if present)
files in the pecans/Mechanisms folder. It passes these files to `mechgen.py` in the pecans/ folder. That
in turn reads those files and generates `chemderiv.pyx`, a Cython source file. `build_pecans.py` then
runs the `setup.py` script in the same directory as itself to convert `chemderiv.pyx` into a `.c` source
file and compile that into a binary executable. All the `chemderiv.*` files will be placed in the pecans/
folder.


## PECANS format

The native PECANS format for the mechanism files is designed to be extremely simple. It is very much modeled
on the Kinetic Pre-processor (KPP) format used in WRF-Chem, but distilled down to its essence.

### Species (.spc) file
List each species used in the mechanism, one per line. Capitalization matters.

### Reaction (.eqn) file
List each reaction, one per line, e.g.:

    NO + O3 = NO2 : ARR2( 1.40e-12, 1310.0, TEMP )

The reaction itself is defined before the `:`. Reactants and products must be separated by _either_ an `=`
or `->`, i.e.

    NO + O3 = NO2

is the same as

    NO + O3 -> NO2

Each reactant and product must be separated by a `+`. A coefficient of 1 may be omitted, as in the previous
examples, while non-unity coefficients must come before the species name to which they apply, e.g.

    1.0 NO + 1.0 O3 = 1.0 NO2

The space between the coefficient and the name is optional. The coefficient must be written so that Python's
`float()` function can understand it; `1`, `1.0`, `1e0`, `1.0e0`, `-1`, etc. are all valid.

The rate constant for the reaction is given after the `:`. It can be one of three forms:

1. A number that `float()` can intepret (see above for the coefficients)
2. A call to a function defined in one of the rates files (next section)
3. A string starting with `j` that indicates a photolysis rate (_not implemented_)

In calls to a defined function, the variables `TEMP` and `CAIR` represent temperature in Kelvin and number
density of air in molecules per cubic centimeter, respectively. E.g.,

    ARR2( 1.40e-12, 1310.0, TEMP )

will call the `ARR2` function with temperature in Kelvin as the third argument.


### Rates (.rate) file

Reaction rates are, unfortunately, one place where it was necessary to require the user to accept a larger
amount of complexity in order to provide sufficient flexibility in defining rate constants. Rate constant
expressions can be defined as Cython `cdef` functions in two places:

1. One of the `.rate` files in pecans/Rates. This is intended for standard rate constants distributed with
the model.
2. A `.rate` file with the same name as the mechanism being compiled, stored in the pecans/Mechanisms folder.
This is intended for user added reactions, although users are welcome to add custom files in the Rates folder
for rate constant expressions which they wish to use across multiple reactions.

Rate expression must be written as valid Cython functions, specifically `cdef` functions (not `def` or `cpdef`,
as `cdef` should be the fastest to call). These look similar to Python function definitions. An example is:

    cdef ko1d(double TEMP, double CAIR):
        cdef double kN
        cdef double k0
    
        kN = 0.78084 * CAIR * 1.8e-11 * exp(107 / TEMP)
        k0 = 0.20946 * CAIR * 3.2e-11 * exp(67 / TEMP)
        return kN + k0

The necessary parts are:
- The keyword `cdef` at the beginning
- The name it will be called by (`ko1d` here)
- The function arguments in parentheses (`TEMP` and `CAIR`), each preceeded by the type `double`.
- A `:` after the closing parenthesis around the arguments
- Any variables used within the function that are not inputs must be declared with the type declaration `cdef double`
- Mathematical expression can be written using normal Python syntax (n.b. that raising to a power is done with `**`,
not `^`). The mathematical functions `exp`, `sqrt`, `log` (natural logarithm), and `log10` (log base 10) are available.
- Return the final value of the rate constant with the `return` keyword.
- For PECANS, temperature in Kelvin and number density of air in molec. / cm^3 must be the variables TEMP and CAIR, 
because the parser looks for those specifically to make sure their position matches in any calls in the `.eqn` file.

As many of these rate expression may be defined in a single `.rate` file as you wish. The parser will raise an exception
if the same name is defined more than once across all `.rate` files read.
