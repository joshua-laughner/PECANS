#!/bin/bash

mydir="$(dirname $0)"
mydir="$(cd $mydir; pwd -P)"
cd "$mydir"

if [[ -f ../bin/activate ]]; then
    echo "Activating python environment at $mydir/../bin/activate"
    source ../bin/activate
else
    echo "Running under system Python"
fi

python -m unittest pecans.transport.tests.simple_transport_tests "$@"
python -m unittest pecans.emissions.tests.ideal_emis_tests "$@"
python -m unittest pecans.chemistry.tests.ideal_chemistry_tests "$@"
# python -m unittest pecans.transport.tests.core_tests "$@"
# python -m unittest pecans.transport.tests.backwards_euler_tests "$@"

echo ""
echo "You may wish to pipe the results into less -S to avoid wrapping long lines"
