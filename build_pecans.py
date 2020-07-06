#!/usr/bin/env python3
import argparse
from distutils.core import run_setup
from glob import glob
import os.path
import sys

from pecans import mechgen

#TODO: add a check that I'm running in a virtual environment?

_mydir = os.path.abspath(os.path.realpath(os.path.dirname(__file__)))
mech_dir = os.path.join(_mydir, 'pecans', 'Mechanisms')

def shell_error(message, exit_code=1):
    print(message, file=sys.stderr)
    _quit(exit_code)

def _make_mechanism_list():
    """
    Make a list of valid mechanisms in the Mechanisms subdirectory of pecans
    :return: A list of valid mechanisms (mechanisms that have both a .eqn and .spc file)
    """
    print(mech_dir)
    mechs_tmp = glob(os.path.join(mech_dir, '*.eqn'))
    mechs_tmp = [os.path.basename(s.replace('.eqn', '')) for s in mechs_tmp]
    mechanisms = [m for m in mechs_tmp if os.path.isfile(os.path.join(mech_dir, m+'.spc'))]
    return mechanisms

def _user_choose_from_list(prompt, options):
    """
    Present the user with a list of available options
    :param options:
    :return:
    """
    print(prompt)
    for opt in options:
        print('  {}: {}'.format(options.index(opt) + 1, opt))

    while True:
        user_ans = input('Enter 1-{} or q to quit: '.format(len(options)))
        if user_ans.lower() == 'q':
            _quit()
        else:
            try:
                ans_index = int(user_ans)-1
            except ValueError:
                print('\n Invalid response.')
            else:
                try:
                    return options[ans_index]
                except IndexError:
                    print('\n Response out of range.')

def _build_mechanism(mechanism, mechanism_style, params):

    species_file = os.path.join(mech_dir, mechanism + '.spc')
    reactions_file = os.path.join(mech_dir, mechanism + '.eqn')
    extra_rate_file = os.path.join(mech_dir, mechanism + '.rate')
    if not os.path.isfile(extra_rate_file):
        extra_rate_file = None

    print('Building {}'.format(mechgen.derivative_file))
    mechgen.generate_mechanism(mechanism_style=mechanism_style, species_file=species_file,
                               reactions_file=reactions_file,
                               additional_rates_file=extra_rate_file, additional_params=params)

    run_setup(os.path.join(_mydir, 'setup.py'), ['build_ext', '--inplace'])

def _quit(exit_code=0):
    if exit_code == 0:
        print('Goodbye.')
    exit(exit_code)

def _get_args():
    parser = argparse.ArgumentParser(description='Build one of the chemical mechanisms for the PECANS model',
                                     epilog='If no mechanism name is provided, a user interactive prompt is given.')
    parser.add_argument('--style', '-s', default='pecans', help='What style the mechanism is. Default is "pecans"')
    parser.add_argument('mechanism', nargs='?', help='The name of the mechanism to use. Must be the base name of and'
                                                     ' .eqn and .spc file in {}'.format(mech_dir))
    parser.add_argument('--params', '-p', metavar="KEY=VALUE", nargs='+', help='Addtional params, the format is a=b.')
    args = parser.parse_args()
    return args


def main():
    print('Running Python from', sys.executable)
    args = _get_args()
    mechanisms = _make_mechanism_list()
    if args.mechanism is None:
        mechanism_to_build = _user_choose_from_list('Choose which mechanism to use:', mechanisms)
    elif args.mechanism not in mechanisms:
        shell_error('"{}" is not a valid mechanism (.spc or .eqn file not found)'.format(args.mechanism))
    else:
        mechanism_to_build = args.mechanism

    _build_mechanism(mechanism_to_build, 'pecans', args.params)
    _quit()

if __name__ == '__main__':
    main()