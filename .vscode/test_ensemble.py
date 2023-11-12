from argparse import ArgumentParser
from functools import partial
from pathlib import Path
from pecans.ensembles import api


_EXAMPLE_DIR = Path(__file__).parent.parent / 'examples'
_TEST_OUTPUT_DIR = Path(__file__).parent.parent / 'test_output'

def debug_ideal(ensemble_mode):
    example_file = _EXAMPLE_DIR / 'one_box_ideal' / 'one_box_ideal.toml'
    ensemble = api.EnsembleRunner(
        base_config_file=str(example_file),
        ensemble_mode=ensemble_mode,
        root_output_dir=str(_TEST_OUTPUT_DIR)
    )
    ensemble.add_ens_var_by_string('CHEMISTRY/mechanism_opts/lifetime_seconds', [3600, 7200, 10800])
    ensemble.add_ens_var_by_string('CHEMISTRY/initial_cond/0/concentration', [2, 20, 200])

    ensemble.run()


def main():
    debug_cases = {
        'ideal_iterations': partial(debug_ideal, 'iterations'),
        'ideal_combinations': partial(debug_ideal, 'combinations'),
    }
    allowed_cases = tuple(debug_cases.keys()) + ('all',)

    p = ArgumentParser(description='Run ensemble test cases in debug mode')
    p.add_argument('test_case', choices=allowed_cases, help='Which case(s) to run')
    clargs = p.parse_args()

    if clargs.test_case == 'all':
        for key, case in debug_cases.items():
            print(f'Running {key}')
            case()
    else:
        print(f'Running {clargs.test_case}')
        case = debug_cases[clargs.test_case]
        case()


if __name__ == '__main__':
    main()
