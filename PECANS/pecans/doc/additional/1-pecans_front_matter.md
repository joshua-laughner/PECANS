# Preface 
## PECANS: the Python Editable Chemical Atmospheric Numeric Solver

The goal of the PECANS multi-box model is to provide a relatively straightforward 
but efficient and flexible idealized atmospheric chemistry modeling framework.
It is not intended to supplant global or regional chemical transport models such
as GEOS-Chem, WRF-Chem, or CMAQ, but instead to offer the capability to carry out
one box to 3D multi-box modeling with idealized (rather than real world) transport.

Keeping the code base clean and easy to follow is a top priority for this model.
This is why we have chosen to develop the code in Python/Cython rather than Fortran
or C, and we will be aggressively documenting the behavior of the model. The code
will also be written with the "code is written once, read many times" philosophy.

Most users should be able to use PECANS in their research without modifying the code
at all, by modifying the mechanism and configuration files. However, this guide will
be written with the expectation that at least some users will want to inspect the code,
either to understand how it works, troubleshoot the behavior of the model, or to
extend the behavior in some manner. Thus, the later chapters of this documentation
will include information about the progammatic structure of the model, as well as
the necessary information for end users.

## Obtaining
The PECANS model can be obtained from the GitHub repository at

<https://github.com/firsttempora/PECANS>

We recommend that you clone the repository, this will allow you to receive updates
most easily. However, you may also download one of the [releases](https://github.com/firsttempora/PECANS/releases).

## Reporting problems or contacting us
The best way to report a problem is via the [Issues](https://github.com/firsttempora/PECANS/issues)
tab of the GitHub repository. This ensures a record of the issue exists and will not be
lost in an inbox somewhere.

If you have a question, it is best if you open an issue with the "question" tag.
This way, there is a record of your question that may help another user and if
it needs to be addressed through an update to the code, it will already be logged.

If you prefer to contact the author(s) directly, reach out to:

[first.tempora@gmail.com](mailto:first.tempora@gmail.com) 

## Contributing

Contributions are welcome! If you would like to submit an improvement, please take the
following steps:

1. Read the style guide and follow the coding style laid out therein.
2. Consider if the improvement keeps the model _straightforward_ and _flexible_ (please contact us if you are uncertain if it does.)
3. [Fork](https://help.github.com/articles/fork-a-repo/) the GitHub repository to your
GitHub account and clone the fork.
4. Make your modifications to the forked repository.
5. Submit a [pull request](https://help.github.com/articles/about-pull-requests/)
