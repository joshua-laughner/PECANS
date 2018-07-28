import netCDF4 as ncdf

from ..emissions import emissions_setup
from ..chemistry import chem_setup
from ..transport import transport_setup
from ..utilities import domain_utilities, io_utils
from ..utilities.config import get_domain_size_from_config


class Domain(object):
    """
    Object that represents the model domain that the simulation occurs on

    :param config: a BetterConfig instance that contains the desired model setup as described in a PECANS config
            file.
    :type config: :class:`~pecans.utilities.config.BetterConfig`

    In a multibox model, the "domain" is the collection of boxes that collectively make up the area that the model is to
    simulate. In PECANS, that is represented by this object. Mainly, this stores information about the size and shape of
    the domain (i.e. how many boxes in each dimension, how large each box it, etc.) and the concentration of all the
    chemical species being tracked. It will also connect to the appropriate functions that calculate the change in
    concentrations due to chemistry, transport, and emissions.

    The Domain object is initialized by calling it with a BetterConfig object from the utilities subpackage. Then to
    advance the model in time, call the `step()` method on the Domain instance, e.g.::

        config = utils.config.load_config_file('pecans_config.cfg')
        model_domain = Domain(config)
        for t in range(100):
            model_domain.step()

    By default, the domain will automatically write an output netCDF file based on the frequency in the configuration
    file. If you want to manually write the model state for any reason, you can call the ``write_output`` method.
    """
    @property
    def seconds_since_model_start(self):
        return self._seconds_since_model_start

    @property
    def species(self):
        return self._chemical_species

    def __init__(self, config):
        self._options = config.section_as_dict('DOMAIN')
        self._config = config
        self._chem_solver, species = chem_setup.setup_chemistry(config)
        self._setup_species(species)
        self._transport_solver, self._get_current_transport = transport_setup.setup_transport(config)
        self._emissions_solver = emissions_setup.setup_emissions(config)
        self._emissions = dict()

        self._seconds_since_model_start = 0

    def _setup_species(self, species):
        """
        Helper method that initializes self._chemical_species based on the configured initial conditions method

        :param species: an iterable of chemical species names
        :type species: iterable of str
        :return: none
        """
        self._chemical_species = dict()
        for specie in species:
            self._chemical_species[specie] = chem_setup.get_initial_conditions(self._config, specie)

    def execute(self, n_seconds_to_run=None):
        """
        Carry out the entire model run

        Once the domain is configured, this will time step the model until the stop time is reached.

        :param n_seconds_to_run: optional, gives the number of seconds that the model should run for. If not given, then
            the run time option in the config is used (which is the most common approach; it is rare that a custom
            number of seconds to run would need to be specified.
        :type n_seconds_to_run: int or float

        :return: none
        """
        if n_seconds_to_run is None:
            n_seconds_to_run = self._options['run_time']

        while self.seconds_since_model_start < n_seconds_to_run:
            self.step()

    def step(self):
        """
        Execute one model time step

        This wil call, in sequence, the configured chemistry solver, transport solver, and emissions solver. It will
        automatically write an output file if the time elapsed since the last output file is longer than the output
        frequency specified in the configuration.

        :return: none
        """
        dt = self._options['dt']
        dx = self._options['dx']
        dy = self._options['dy']
        dz = self._options['dz']
        domain_size = get_domain_size_from_config(self._config)

        # Start by creating the delta array with the chemistry, since that will automatically set up delta as a
        # dictionary with the same keys as self._chemical_species
        if self._config.get('CHEMISTRY', 'do_chemistry'):
            self._chemical_species = self._chem_solver(dt, TEMP=None, CAIR=None, **self._chemical_species)

        # Now we need to handle emissions and transport
        if self._config.get('TRANSPORT', 'do_transport'):
            for name, concentration in self._chemical_species.items():
                u_x, u_y, u_z, D_x, D_y, D_z = self._get_current_transport(self._seconds_since_model_start)
                # This is one way of handling the transport code, which assumes that dy, dz are set to None if the model
                # doesn't have that dimension. A better way would be to rework the transport code to be smarter about
                # how it checks the input, and either set the unnecessary values to None inside it, or just ignore them
                dy_tmp = dy if len(domain_size) >= 2 else None
                dz_tmp = dz if len(domain_size) >= 3 else None
                self._chemical_species[name] = self._transport_solver(concentration, dt, dx=dx, dy=dy_tmp, dz=dz_tmp,
                                                                      u_x=u_x, u_y=u_y, u_z=u_z, D_x=D_x, D_y=D_y, D_z=D_z,
                                                                      domain_size=domain_size)

        if self._config.get('EMISSIONS', 'do_emissions'):
            self._chemical_species, self._emissions = self._emissions_solver(config=self._config,
                                                                             seconds_since_model_start=self.seconds_since_model_start,
                                                                             **self._chemical_species)

        output_frequency = self._config.get('OUTPUT', 'output_frequency')
        self._seconds_since_model_start += dt
        if output_frequency != 0 and self._seconds_since_model_start % output_frequency < dt:
            self.write_output()

        return self.seconds_since_model_start

    def write_output(self, output_file_name=None):
        """
        Write an output netCDF file representing the instantaneous current model state

        :param output_file_name: optional, allows you to specify the desired output file name. If not given, it will
            default to

            "pecans_output_DDDdHHhMMmSSs.nc"

            where DDD, HH, MM, and SS are the days, minutes and seconds since the beginning of the model run.
        :type output_file_name: str

        :return: none
        """
        def create_dim_and_coord(dataset, name, coordinates):
            dataset.createDimension(name, len(coordinates))
            coord_var = dataset.createVariable(name, io_utils.data_type, (name,))
            coord_var[:] = coordinates

        if output_file_name is None:
            seconds = self._seconds_since_model_start
            days = self._seconds_since_model_start // (60 * 60 * 24)
            seconds -= days * (60 * 60 * 24)
            hours = seconds // 3600
            seconds -= hours * 3600
            minutes = seconds // 60
            seconds -= minutes * 60

            output_file_name = 'pecans_output_{:03}d{:02}h{:02}m{:02}s.nc'.format(days, hours, minutes, seconds)

        print('Writing {}'.format(output_file_name))
        with ncdf.Dataset(output_file_name, mode='w', clobber=True, format='NETCDF4') as ncdat:
            x_coord, y_coord, z_coord = domain_utilities.compute_coordinates_from_config(config=self._config)

            create_dim_and_coord(ncdat, 'x', x_coord)
            dims = ['x']

            if y_coord is not None:
                create_dim_and_coord(ncdat, 'y', y_coord)
                dims.append('y')

            if z_coord is not None:
                create_dim_and_coord(ncdat, 'z', z_coord)
                dims.append('z')

            for name, conc in self._chemical_species.items():
                this_var = ncdat.createVariable(name, io_utils.data_type, tuple(dims))
                this_var[:] = conc

            for name, emis in self._emissions.items():
                this_var = ncdat.createVariable('E_' + name, io_utils.data_type, tuple(dims))
                this_var[:] = emis
