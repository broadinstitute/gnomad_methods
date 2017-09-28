import unittest

from utils import *

hc = None
verbose = False
test_resources_dir = 'tests/resources'


def setUpModule():
    global hc
    global verbose
    verbose = '-v' in sys.argv
    hc = HailContext(log='/dev/null', master='local[1]')


def tearDownModule():
    global hc
    hc.stop()
    hc = None


class FilteringTests(unittest.TestCase):

    @staticmethod
    def create_filter_test_vds():
        """

        :return: VDS with some filters
        :rtype: VariantDataset
        """
        rows = [
            # Bi-allelic expected behavior
            {'v': Variant.parse('1:10000:A:T'),   'InbreedingCoeff': None, 'AS_FilterStatus': [[]],              'expected_filters': [],                        'expected_after_split': [[]]},
            {'v': Variant.parse('1:10001:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[]],              'expected_filters': [],                        'expected_after_split': [[]]},
            {'v': Variant.parse('1:10002:A:T'),   'InbreedingCoeff': -0.5, 'AS_FilterStatus': [[]],              'expected_filters': ['InbreedingCoeff'],       'expected_after_split': [['InbreedingCoeff']]},
            {'v': Variant.parse('1:10003:A:T'),   'InbreedingCoeff': -0.5, 'AS_FilterStatus': [['RF']],          'expected_filters': ['InbreedingCoeff', 'RF'], 'expected_after_split': [['InbreedingCoeff', 'RF']]},
            {'v': Variant.parse('1:10004:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF']],          'expected_filters': ['RF'],                    'expected_after_split': [['RF']]},
            {'v': Variant.parse('1:10005:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF', 'AC0']],   'expected_filters': ['RF', 'AC0'],             'expected_after_split': [['RF', 'AC0']]},

            # Multi-allelic expected behavior
            {'v': Variant.parse('2:10000:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[], []],          'expected_filters': [],                               'expected_after_split': [[], []]},
            {'v': Variant.parse('2:10001:A:T,C'), 'InbreedingCoeff': -0.5, 'AS_FilterStatus': [[], []],          'expected_filters': ['InbreedingCoeff'],              'expected_after_split': [['InbreedingCoeff'], ['InbreedingCoeff']]},
            {'v': Variant.parse('2:10002:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF'], []],      'expected_filters': [],                               'expected_after_split': [['RF'], []]},
            {'v': Variant.parse('2:10003:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF'], ['RF']],  'expected_filters': ['RF'],                           'expected_after_split': [['RF'], ['RF']]},
            {'v': Variant.parse('2:10004:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF'], ['AC0']], 'expected_filters': ['RF', 'AC0'],                    'expected_after_split': [['RF'], ['AC0']]},
            {'v': Variant.parse('2:10005:A:T,C'), 'InbreedingCoeff': -0.5, 'AS_FilterStatus': [['RF'], []],      'expected_filters': ['InbreedingCoeff'],              'expected_after_split': [['InbreedingCoeff', 'RF'], ['InbreedingCoeff']]},
            {'v': Variant.parse('2:10006:A:T,C'), 'InbreedingCoeff': -0.5, 'AS_FilterStatus': [['RF'], ['AC0']], 'expected_filters': ['InbreedingCoeff', 'RF', 'AC0'], 'expected_after_split': [['InbreedingCoeff', 'RF'], ['InbreedingCoeff', 'AC0']]},

            # Unexpected behavior
            {'v': Variant.parse('9:10000:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': None,              'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10001:A:T'),   'InbreedingCoeff': None, 'AS_FilterStatus': None,              'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10002:A:T'),   'InbreedingCoeff': -0.5, 'AS_FilterStatus': None,              'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10003:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [None],            'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10004:A:T'),   'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[None]],          'expected_filters': None,                      'expected_after_split': None},
            {'v': Variant.parse('9:10005:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[], None],        'expected_filters': None,                      'expected_after_split': [[], None]},
            {'v': Variant.parse('9:10006:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [[], [None]],      'expected_filters': None,                      'expected_after_split': [[], None]},
            {'v': Variant.parse('9:10007:A:T,C'), 'InbreedingCoeff': 0.0,  'AS_FilterStatus': [['RF'], [None]],  'expected_filters': None,                      'expected_after_split': [['RF'], None]},
        ]
        schema = ['v', 'InbreedingCoeff', 'AS_FilterStatus', 'expected_filters', 'expected_after_split']
        types = [TVariant(), TDouble(), TArray(TSet(TString())), TSet(TString()), TArray(TSet(TString()))]
        return VariantDataset.from_table(KeyTable.from_py(hc, rows, TStruct(schema, types), key_names=['v']))

    @classmethod
    def setUpClass(cls):
        cls.vds = cls.create_filter_test_vds()
        if verbose: cls.vds.variants_table().show(50)

    def test_allele_filtering(self):
        site_filters = {
            'InbreedingCoeff': 'isDefined(va.InbreedingCoeff) && va.InbreedingCoeff < -0.3'
        }

        result_vds = set_site_filters(self.vds, site_filters, 'va.AS_FilterStatus')
        if verbose: result_vds.variants_table().show(50)
        result = result_vds.query_variants('variants.map(v => (isMissing(va.filters) && isMissing(va.expected_filters)) || va.filters == va.expected_filters).counter()')
        self.assertEqual(result[True], sum(result.values()))

        split_vds = result_vds.split_multi().annotate_variants_expr(index_into_arrays(['va.AS_FilterStatus', 'va.expected_after_split']))
        result_split_vds = set_site_filters(split_vds, site_filters, 'va.AS_FilterStatus')
        if verbose: result_split_vds.variants_table().show(50)
        result = result_split_vds.query_variants('variants.map(v => (isMissing(va.filters) && isMissing(va.expected_filters)) || va.filters == va.expected_after_split).counter()')
        self.assertEqual(result[True], sum(result.values()))


class KeyTableTests(unittest.TestCase):

    @staticmethod
    def create_frequency_kt():
        """
        KeyTable with some frequency data

        :return: keytable with frequency data
        :rtype: KeyTable
        """
        rows = [
            {'v': Variant.parse('1:10000:A:T'), 'AC_NFE': 1,  'AC_AFR': 8,   'Hom_NFE': 0, 'Hom_AFR': 0},
            {'v': Variant.parse('1:10001:A:T'), 'AC_NFE': 10, 'AC_AFR': 100, 'Hom_NFE': 1, 'Hom_AFR': 10},
        ]
        schema = ['v', 'AC_NFE', 'AC_AFR', 'Hom_NFE', 'Hom_AFR']
        types = [TVariant(), TInt(), TInt(), TInt(), TInt()]
        return KeyTable.from_py(hc, rows, TStruct(schema, types), key_names=['v'])

    @classmethod
    def setUpClass(cls):
        cls.kt = cls.create_frequency_kt()
        if verbose: cls.kt.show(50)

    def test_melt_kt(self):
        melted_kt = melt_kt(self.kt, columns_to_melt=['AC_NFE', 'AC_AFR', 'Hom_NFE', 'Hom_AFR'])
        self.assertEqual(melted_kt.count(), 8)
        self.assertEqual(sorted(melted_kt.columns), sorted(['v', 'value', 'variable']))
        self.assertEqual(melted_kt.query('variable.counter()'), {'AC_NFE': 2, 'AC_AFR': 2, 'Hom_NFE': 2, 'Hom_AFR': 2})
        self.assertEqual(melted_kt.filter('v == Variant("1:10000:A:T") && variable == "AC_NFE"').query('value.collect()'), [1])
        if verbose: melted_kt.show(50)

    def test_melt_grouped_kt(self):
        grouped_melted_kt = melt_kt_grouped(self.kt,
                                            columns_to_melt={'NFE': ['AC_NFE', 'Hom_NFE'], 'AFR': ['AC_AFR', 'Hom_AFR']},
                                            value_column_names=['AC', 'Hom'],
                                            key_column_name='pop')
        self.assertEqual(grouped_melted_kt.count(), 4)
        self.assertEqual(sorted(grouped_melted_kt.columns), sorted(['v', 'pop', 'AC', 'Hom']))
        self.assertEqual(grouped_melted_kt.query('pop.counter()'), {'NFE': 2, 'AFR': 2})
        if verbose: grouped_melted_kt.show(50)


class VEPTests(unittest.TestCase):

    @staticmethod
    def generate_stress_test_vds():
        rows = [
            {'v': Variant.parse('1:69270:A:G'),       'test_transcript': 'ENST00000335137', 'expected_lof': '',   'expected_lof_filter': '',                    'expected_lof_flags': ''},
            {'v': Variant.parse('1:69869:T:A'),       'test_transcript': 'ENST00000335137', 'expected_lof': 'HC', 'expected_lof_filter': '',                    'expected_lof_flags': 'SINGLE_EXON'},
            {'v': Variant.parse('1:739130:T:TA'),     'test_transcript': 'ENST00000599533', 'expected_lof': 'LC', 'expected_lof_filter': 'NON_CAN_SPLICE_SURR', 'expected_lof_flags': 'PHYLOCSF_WEAK'},
            {'v': Variant.parse('1:861301:G:A'),      'test_transcript': 'ENST00000342066', 'expected_lof': 'HC', 'expected_lof_filter': '',                    'expected_lof_flags': ''},
            {'v': Variant.parse('1:896934:C:T'),      'test_transcript': 'ENST00000338591', 'expected_lof': 'LC', 'expected_lof_filter': 'NON_CAN_SPLICE',      'expected_lof_flags': ''},
            {'v': Variant.parse('1:900341:A:C'),      'test_transcript': 'ENST00000338591', 'expected_lof': 'HC', 'expected_lof_filter': '',                    'expected_lof_flags': 'NAGNAG_SITE'},
            {'v': Variant.parse('1:915034:C:T'),      'test_transcript': 'ENST00000433179', 'expected_lof': 'HC', 'expected_lof_filter': '',                    'expected_lof_flags': 'PHYLOCSF_WEAK'},
            {'v': Variant.parse('1:1018307:G:T'),     'test_transcript': 'ENST00000434641', 'expected_lof': 'LC', 'expected_lof_filter': 'END_TRUNC',           'expected_lof_flags': 'PHYLOCSF_UNLIKELY_ORF'},
            {'v': Variant.parse('1:1226966:GC:G'),    'test_transcript': 'ENST00000379116', 'expected_lof': 'LC', 'expected_lof_filter': 'END_TRUNC',           'expected_lof_flags': ''},  # should fail 50_BP_RULE
            {'v': Variant.parse('1:1341007:G:T'),     'test_transcript': 'ENST00000482352', 'expected_lof': 'LC', 'expected_lof_filter': 'END_TRUNC',           'expected_lof_flags': ''},  # should fail 50bp filter
            {'v': Variant.parse('1:1653047:G:C'),     'test_transcript': 'ENST00000378638', 'expected_lof': 'LC', 'expected_lof_filter': 'ANC_ALLELE',          'expected_lof_flags': ''},
            {'v': Variant.parse('1:3647537:C:T'),     'test_transcript': 'ENST00000378280', 'expected_lof': 'HC', 'expected_lof_filter': '',                    'expected_lof_flags': 'PHYLOCSF_UNLIKELY_ORF'},  # Should pass 50 bp filter (in 2nd to last exon, but past 50 bp)
            {'v': Variant.parse('1:9656068:G:GGTGT'), 'test_transcript': 'ENST00000340305', 'expected_lof': 'LC', 'expected_lof_filter': 'EXON_INTRON_UNDEF',   'expected_lof_flags': ''},
        ]

    @classmethod
    def setUpClass(cls):
        cls.vds = hc.read('{}/{}'.format(test_resources_dir, 'vep_test.vds'))

    def process_consequences_test(self):
        proc_vds = process_consequences(self.vds)


if __name__ == '__main__':
    unittest.main()