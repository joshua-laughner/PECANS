import os

from pecans.utilities.config import load_config_file
from pecans.core import Domain

import pdb


def main():
    config_file_name = os.path.join(os.path.dirname(__file__), 'pecans_config.cfg')
    config = load_config_file(config_file_name)
    model_domain = Domain(config)
    curr_t = model_domain.seconds_since_model_start
    while curr_t < config.get('DOMAIN', 'run_time'):
        curr_t = model_domain.step()


if __name__ == '__main__':
    main()