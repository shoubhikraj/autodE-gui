# Contains the config class which is instantiated once at the beginning
# of the

class _AdeGuiConfig:
    ade_n_cores = 8  # Number of cores in autodE
    ade_job_typ = 'Reaction Profile'

    # TODO: remove later <= this is not needed as we are in GUI
    def __setattr__(self, key, value):
        """ Check if the values are allowed """
        if not hasattr(self, key):
            raise Exception(f'Cannot set {key}, it is not present in autodE-GUI')

        allowed_job_typs = ['Reaction Profile','Reaction Path: NEB','Reaction Path: CI-NEB',
                            'Reaction Path: Dimer', '1D PES', '2D PES']
        if key == 'ade_job_typ':
            if value not in allowed_job_typs:
                raise Exception(f'{key} is not an allowed job type')


Config = _AdeGuiConfig()
