import netCDF4 as ncdf

from ..chemistry import chem_setup, emissions_setup
from ..transport import transport_setup
from ..utilities.config import ConfigurationError, get_domain_size_from_config
from ..utilities import domain_utilities, io_utils

import pdb


class Domain(object):
    @property
    def seconds_since_model_start(self):
        return self._seconds_since_model_start

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
        self._chemical_species = dict()
        for specie in species:
            self._chemical_species[specie] = chem_setup.get_initial_conditions(self._config, specie)

    def step(self):
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
            self._emissions = dict() # clear previous timestep's emissions
            for name, concentration in self._chemical_species.items():
                #pdb.set_trace()
                self._chemical_species[name], self._emissions[name] = self._emissions_solver(self._config, concentration, name,
                                                                                             self.seconds_since_model_start)

        self._seconds_since_model_start += dt
        if self._seconds_since_model_start % self._config.get('OUTPUT', 'output_frequency') < dt:
            self.write_output()

        return self.seconds_since_model_start

    def write_output(self):
        def create_dim_and_coord(dataset, name, coordinates):
            dataset.createDimension(name, len(coordinates))
            coord_var = dataset.createVariable(name, io_utils.data_type, (name,))
            coord_var[:] = coordinates

        seconds = self._seconds_since_model_start
        days = self._seconds_since_model_start // (60 * 60 * 24)
        seconds -= days * (60 * 60 * 24)
        hours = seconds // 3600
        seconds -= hours * 3600
        minutes = seconds // 60
        seconds -= minutes * 60

        output_file_name = 'pecans_output_{:03}d{:02}h{:02}m{:02}s.nc'.format(days, hours, minutes, seconds)
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
