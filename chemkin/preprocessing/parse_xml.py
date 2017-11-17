from enum import Enum
import xml.etree.ElementTree as ET
import numpy as np

from chemkin.chemkin_errors import ChemKinError
from chemkin.reaction.reaction_coefficients import ArrheniusCoefficient, \
    ConstantCoefficient, ModifiedArrheniusCoefficient
from chemkin.thermodynamics.thermo import Thermo


class RxnType(Enum):
    Elementary = 1


class XmlParser():
    """Core method load() produces list of RxnData from XML file contents.

    Notes:
        ChemKinError raised when invalid values are encountered in the XML file.
    """

    def __init__ (self, path):
        """ Ensures path contains .xml file extension. """
        if path[-4:] != '.xml':
            path += '.xml'
        self.path = path

    def load (self):
        """ Parses XML file contents to create list of RxnData objects
        representing the reactions in the file.
        """
        tree = ET.parse(self.path)
        root = tree.getroot()
        species = []
        for species_i in root.find('phase'):
            new_species = species_i.text.strip().split()
            species.extend(new_species)

        results = []
        for rxn in root.find('reactionData').findall('reaction'):
            rxn_data = self.__extract_data_from_reaction_element(rxn)
            results.append(rxn_data)
        return species, results

    def __extract_data_from_reaction_element (self, rxn):
        """Returns RxnData object containing data from <reaction> XML element.
        Raises ChemKinError for invalid attribute/element values.
        """
        result = RxnData()
        result.rxn_id = rxn.get('id')

        # reversible
        reversible = rxn.get('reversible').lower().strip()
        if reversible in ['no', 'n', 'false', 'f']:
            result.reversible = False
        elif reversible in ['yes', 'y', 'true', 't']:
            result.reversible = True
        else:
            raise ChemKinError(
                  'XmlParser.load()',
                  'Invalid reversible attribute in reaction {}'.format(
                        result.rxn_id))

        # type
        rnx_type = rxn.get('type').lower().strip()
        if rnx_type == 'elementary':
            result.type = RxnType.Elementary
        else:
            raise ChemKinError(
                  'XmlParser.load()',
                  'Invalid type attribute in reaction {}'.format(
                        result.rxn_id))

        # rate_coeff
        rate_coeff = rxn.find('rateCoeff')
        if rate_coeff is None:
            raise ChemKinError('XmlParser.load()',
                               'No <rateCoeff> element found in one of the '
                               'reactions.')

        if rate_coeff.find('Arrhenius') is not None:
            arrhenius = rate_coeff.find('Arrhenius')
            A = float(arrhenius.find('A').text.strip())
            if A < 0:
                raise ChemKinError('XmlParser.load()',
                                   'A coeff < 0 in reaction with ' \
                                   'id = {}'.format(result.rxn_id))
            E = float(arrhenius.find('E').text.strip())
            result.rate_coeff = [A, E]
        elif rate_coeff.find('modifiedArrhenius') is not None:
            mod_arrhenius = rate_coeff.find('modifiedArrhenius')
            A = float(mod_arrhenius.find('A').text.strip())
            if A < 0:
                raise ChemKinError(
                      'A coeff < 0 in reaction with id = {}'.format(
                            result.rxn_id))
            b = float(mod_arrhenius.find('b').text.strip())
            E = float(mod_arrhenius.find('E').text.strip())
            result.rate_coeff = [A, b, E]
        elif rate_coeff.find('Constant') is not None:
            const = rate_coeff.find('Constant')
            result.rate_coeff = float(const.find('k').text.strip())
        else:
            raise ChemKinError('XmlParser.load()',
                               'No recognized child of <rateCoeff> found '
                               'from which to parse coefficients.')

        # reactants / products
        result.reactants = self.__map_conc_to_species(rxn.find('reactants'))
        result.products = self.__map_conc_to_species(rxn.find('products'))

        return result

    @staticmethod
    def __map_conc_to_species (tag):
        """ Creates dict mapping the species to concentration where tag is
        either a <reactants> or <products> XML tag.
        """
        result = {}
        for item in tag.text.strip().split():
            species, conc = item.strip().split(':')
            conc = int(conc)
            species = species.upper()
            result[species] = conc
        return result

    def parsed_data_list (self, Ti):
        species, rxn_data_list = self.load()
        n_species = len(species)

        species_idx_dict = {}  # build the dictionary of key = species_name,
        # value = species_index
        for i, s in enumerate(species):
            species_idx_dict[s] = i

        parsed_data_dic_list = []  # list of dicts, one for each 1 set of
        # rxns (1 xml file) under one temperature
        for T in Ti:
            sys_vi_p = []  # list of reactant Stoichiometric coefficients in
            # each rxn
            sys_vi_dp = []  # list of product Stoichiometric coefficients in
            # each rxn
            ki = []  # list of reation rate coefficients in each rxn
            is_reversible = None  # indicator of the system of reactions
            # being irreversible/reversible

            for rxn_data in rxn_data_list:  # 1 rxn per rxn_data
                if rxn_data.type != RxnType.Elementary:
                    raise ChemKinError('XmlParser.parsed_data_list(Ti)',
                                       'Non-elementary reactions cannot be '
                                       'parsed now.')

                if is_reversible == None:
                    is_reversible = rxn_data.reversible  # set the indicator
                    # of the system of reactions to be irreversible/reversible

                if rxn_data.reversible != is_reversible:  # the system of
                    # reactions in the XML file must be all
                    # irreversible/reversible
                    raise ChemKinError('XmlParser.parsed_data_list(Ti)',
                                       'The system of reactions in the XML '
                                       'file {} are inconsistent in '
                                       'reversibility.'.format(
                                           self.path))

                rxn_id = rxn_data.rxn_id  # save id

                rxn_vi_p = np.zeros((
                                    n_species,))  # save the Stoichiometric
                # coefficients of the reactants in this rxn
                for s, vi in rxn_data.reactants.items():
                    idx = species_idx_dict[s]  # get index of the specii
                    rxn_vi_p[idx] = vi
                sys_vi_p.append(list(rxn_vi_p))

                rxn_vi_dp = np.zeros((
                                     n_species,))  # save the Stoichiometric
                # coefficients of the products in this rxn
                for s, vi in rxn_data.products.items():
                    idx = species_idx_dict[s]  # get index of the specii
                    rxn_vi_dp[idx] = vi
                sys_vi_dp.append(list(rxn_vi_dp))

                coef_params = rxn_data.rate_coeff
                if isinstance(coef_params, list):
                    if len(coef_params) == 3:  # modified arrhenius coef
                        A = coef_params[0]
                        b = coef_params[1]
                        E = coef_params[2]
                        ki.append(
                            ModifiedArrheniusCoefficient(A, b, E, T).get_coef())
                    else:  # arrhenius coef
                        A = coef_params[0]
                        E = coef_params[1]
                        ki.append(ArrheniusCoefficient(A, E, T).get_coef())
                else:  # const coef
                    ki.append(ConstantCoefficient(coef_params).get_coef())

                parsed_data_dic = {}
                parsed_data_dic['species'] = species
                parsed_data_dic['ki'] = ki
                parsed_data_dic['sys_vi_p'] = sys_vi_p
                parsed_data_dic['sys_vi_dp'] = sys_vi_dp
                parsed_data_dic['is_reversible'] = is_reversible
                parsed_data_dic['T'] = T

            if is_reversible == True:
                try:
                    b_ki = Thermo(species, T, ki, sys_vi_p,
                                  sys_vi_dp).get_backward_coefs()
                    parsed_data_dic['b_ki'] = b_ki
                except ChemKinError as err:
                    print(str(err))
                    parsed_data_dic['b_ki'] = 'Not Defined'

            parsed_data_dic_list.append(parsed_data_dic)

        return parsed_data_dic_list


class RxnData():
    """ Container for individual reaction data.

    Attributes
    ----------
    rxn_id : str
        id attribute of <reaction> element.
    reversible : bool
        True if reaction is reversible; False if irreversible.
    reactants : Dict[str, int]
        Mapping of species to stoichiometric coefficients for reactants.
        Example: {'H2':1, 'O':1}
    products : Dict[str, int]
        Mapping of species to stoichiometric coefficients for products.
    rate_coeff : List[float] or float
        Reaction rate coefficients contained depend on the type of rate
        coefficients in the XML file, dictated by the child of the <rateCoeff>
        element:
            Arrhenius: [A, E]
            modifiedArrhenius: [A, b, E]
            Constant: k
    type : RxnType
        Enum value for reaction type.
    """

    def __init__ (self, rxn_id=None, reversible=None, reactants=None,
                  products=None, rate_coeff=None, type=None):
        self.rxn_id = rxn_id
        self.reversible = reversible
        self.reactants = reactants
        self.products = products
        self.rate_coeff = rate_coeff
        self.type = type

    def equation (self):
        """ Returns equation representation of reactants and products.

        Species are listed in alphabetical order on both reactant and product
        side.
        """
        reactant_side = self.__build_equation_side(self.reactants)
        product_side = self.__build_equation_side(self.products)
        return '{0} [=] {1}'.format(reactant_side, product_side)

    @staticmethod
    def __build_equation_side (species_dict):
        result = ''
        species_count = -1
        for species, conc in sorted(species_dict.items()):
            species_count += 1

            # Handle first species differently than rest, then proceed to
            # next species.
            if species_count == 0:
                if abs(conc) == 1:
                    conc_part = ''
                else:
                    conc_part = conc
                result += '{0}{1}'.format(conc_part, species)
                continue

            # Non-first species.
            if conc < 0:
                operator = ' - '
            else:
                operator = ' + '

            if abs(conc) == 1:
                conc_part = ''
            else:
                conc_part = abs(conc)

            result += '{0}{1}{2}'.format(operator, conc_part, species)

        return result
