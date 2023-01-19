    def _init_phase_details(self):

        # Read in allowable phases
        phase_info, info_files = phases.get_phase_info()

        def new_phase_dict():
            new_dict = OrderedDict()
            new_dict['pure'] = None
            new_dict['solution'] = None

            new_dict['active_pure'] = None
            new_dict['active_solution'] = None
            return new_dict

        # Setup phase_details
        phase_details = OrderedDict()
        phase_details['database'] = self.database
        phase_details['filenames'] = new_phase_dict()
        phase_details['info'] = new_phase_dict()

        # load allowable phases into phase_details
        phase_details['filenames']['pure'] = info_files['pure']
        phase_details['info']['pure'] = phase_info['pure']

        phase_details['filenames']['solution'] = info_files['solution']
        phase_details['info']['solution'] = phase_info['solution']

        self._phase_details = phase_details

    def _init_active_phases(self, fixH2O=True, liq_mod=None):
        self._liq_mod = liq_mod

        # phase_details = self.phase_details
        # self._validate_active_phases()

        pure_propsDB, pure_objDB = (
            self._init_database('pure', fixH2O=fixH2O))

        solution_propsDB, solution_objDB = (
            self._init_database('solution', liq_mod=liq_mod))

        self._phase_props = {'pure':pure_propsDB,
                             'solution':solution_propsDB}

        phase_obj = OrderedDict()

        for ikey, ival in pure_objDB.items():
            phase_obj[ikey] = ival

        for ikey, ival in solution_objDB.items():
            phase_obj[ikey] = ival


        self._phase_obj = phase_obj
        # self._phase_obj = {'pure':pure_objDB,
        #                    'solution':solution_objDB}

    def _validate_active_phases(self):
        phase_details = self.phase_details

        def load_validate_store_info(basename, phase_type):
            # load
            phase_info = phase_details['info'][phase_type]

            active_phases, active_filename = self._read_database_phase_file(basename)

            # validate
            abbrev_valid = active_phases['Abbrev'].isin(
                phase_info['Abbrev'])

            err_msg = (
                'The {phase_type} phase library defined in '
                '{filename} contains some invalid phase '
                'abbreviations, shown below: \n\n'
                '{invalid_phases}\n\n'
                'Check that the abbreviations conform to the '
                'list given in: "{phase_info_file}"')

            invalid_phases = str(active_phases[~abbrev_valid])
            phase_info_file = phase_details['filenames'][phase_type]

            assert abbrev_valid.all(), (
                err_msg.format(phase_type=phase_type,
                               filename=active_filename,
                               invalid_phases=invalid_phases,
                               phase_info_file=phase_info_file)
            )

            # store
            phase_details['filenames']['active_'+phase_type] = active_filename
            phase_details['info']['active_'+phase_type] = active_phases
            pass

        load_validate_store_info('PurePhases.csv', 'pure')
        load_validate_store_info('SolutionPhases.csv', 'solution')


    @property
    def phase_details(self):
        """
        Dictionary containing phase information.

        Phase info for global phase list and active phases for this database.
        Dictionary includes:
            'database': str ID of chosen thermodynamic database
            'filenames': dict of filenames where raw phase info is stored
            'info': dict of phase info tables
                'pure': global pure phase info table
                'solution': global solution phase info table
                'active_pure': active pure phase info table for this database
                'active_solution': active solution phase info table for this database
        """
        return self._phase_details


    def _init_phase_attributes(self):
        phase_obj = self._phase_obj
        phase_props = self._phase_props
        phase_details = self._phase_details

        def create_attributes(name=None,
                              phase_type=None,
                              props=None,
                              obj=None,
                              parents=None,
                              members=None):

            attributes = OrderedDict()
            attributes['name'] = name
            attributes['phase_type'] = phase_type
            attributes['parents'] = parents
            attributes['members'] = members
            attributes['props'] = props
            attributes['obj'] = obj
            return attributes

        phase_attributes = OrderedDict()

        for phase_type in phase_obj:
            info = phase_details['info'][phase_type]
            # active_info = phase_details['info']['active_'+phase_type]

            for abbrev in phase_obj[phase_type]:
                mask_i = info['Abbrev']==abbrev
                info_i = info.loc[mask_i]

                attributes = create_attributes(
                    name=info_i['Name'],
                    phase_type=phase_type,
                    obj=phase_obj[phase_type][abbrev],
                    props=phase_props[phase_type][abbrev])

                phase_attributes[abbrev] = attributes

        self._phase_attributes = phase_attributes


        """
        Dictionary of phase attributes, stored by phase abbreviations.

        For each phase, there is a dictionary containing:
            'name':
            'phase_type':
            'obj':
            'props':
                'abbrev':
                'phase_name':
                'class_name':
                'endmember_name':
                'endmember_id':
                'formula':
                'atom_num':
                'mol_wt':
                'elemental_entropy':
                'elemental_comp':
                'mol_oxide_comp':
        """


    def _init_phase_database(self, phase_type, fixH2O=True, liq_mod=None):
        objDB = OrderedDict()
        phase_obj = None
        for indx, phase_info in active_info.iterrows():
            abbrev = phase_info['Abbrev']
            classnm = phase_info['ClassName']
            phase_obj = self._init_phase_obj(
                classnm+self.database, abbrev,
                phase_type)
            try:

                propsDB[abbrev] = phase_obj.props
                objDB[abbrev] = phase_obj

            except:
                print('{classnm} is not a valid ClassName for '
                    'the {database} database.'.\
                    format(classnm=classnm, database=self.database))

        propsDB, objDB = self._load_special_phases(phase_type, fixH2O, liq_mod,
                                                  propsDB, objDB)
        return propsDB, objDB
