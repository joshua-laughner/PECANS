import os
from pecans.utilities.config import load_config_file
from pecans.main import Domain


def main():
    config_file_name = os.path.join(os.path.dirname(__file__), 'pecans_config.toml')
    config = load_config_file(config_file_name)
    model_domain = Domain(config)
    model_domain.execute()


if __name__ == '__main__':
    main()
