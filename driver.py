from chemkin import pckg_xml_path
from chemkin.preprocessing.parse_xml import XmlParser
from chemkin.viz import summary
from chemkin.reaction.elementary_rxn import ReversibleElementaryRxn
from chemkin.reaction.elementary_rxn import IrreversibleElementaryRxn


# Ti = [750, 1500, 2500]
# xi = [2.0, 1.0, 0.5, 1.0, 1.0]  # specie concentrations for 'rxns_hw5.xml'
# xml_parser = XmlParser(pckg_xml_path('rxns_hw5'))

Ti = [1500]
xi = [2., 1., .5, 1., 1., 1., .5, 1.] # specie concentrations 'rxns_reversible.xml'
xml_parser = XmlParser(pckg_xml_path('rxns_reversible'))

parsed_data_list = xml_parser.parsed_data_list(Ti)
# summary.print_reaction_rate(parsed_data_list, xi)

test_flag = 0 # reaction rates can be printed

for parsed_data in parsed_data_list:

    species = parsed_data['species']
    ki = parsed_data['ki']
    sys_vi_p = parsed_data['sys_vi_p']
    sys_vi_dp = parsed_data['sys_vi_dp']
    is_reversible = parsed_data['is_reversible']
    T = parsed_data['T']

    if is_reversible == False:
        progress_rates = IrreversibleElementaryRxn(ki, xi, sys_vi_p, sys_vi_dp).progress_rate()

    else:
        b_ki = parsed_data['b_ki']
        if str(b_ki) == 'Not Defined':
            test_flag = 1  # reaction rates cannot be printed because T is not in some specie's temperature range

            print('------At Temperature', T, 'K------')
            print(
            'Backward reaction coefficients not defined: T={} is not in some specie\'s temperature range.'.format(T))
            print('--------------------------------')
            continue

        progress_rates = ReversibleElementaryRxn(ki, b_ki, xi, sys_vi_p, sys_vi_dp).progress_rate()

        print(progress_rates)
        print(species)
        print(ki)
        print(b_ki)
        print(sys_vi_p)
        print(sys_vi_p)

        ReversibleElementaryRxn([10, 10], [10, 10], [1.0, 2.0, 1.0], [[1.0, 2.0, 0.0], [2.0, 0.0, 2.0]],
                                  [[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]]).progress_rate()