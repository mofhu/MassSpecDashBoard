import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

class Search_Result():
    def __init__(self, filename='', identifier='',
                 search_workflow='ID_PD'):
        """use resultname.filename, resultname.identifier, etc.

        could extend more using a method for i/o.
        also could save in class.
        """
        self.filename = filename  # implies where to read
        self.identifier = identifier
        self.search_workflow = search_workflow  # ID/Quan_search engine (top-down hierarchy)

    def read_result(self):
        # output: self.raw_result, self.result
        # a dict {protein: protein_result, peptide: peptide_result, PSM: PSM_result}
        self.raw_result = {}
        if self.search_workflow == 'ID_PD':
            suffix = {'PSM': self.filename + '_PSMs.txt',
                      'peptide': self.filename + '_PeptideGroups.txt',
                      'protein': self.filename + '_Proteins.txt'}
            for key in suffix:
                self.raw_result[key] = pd.read_csv(suffix[key], sep = '\t')
                self.raw_result[key]['Identifier'] = self.identifier

        # raw result to final result
        self.result = {}
        if self.search_workflow == 'ID_PD':
            # select columns
            # protein: no change, only filter columns
            protein_cols = ['Identifier', 'Accession', 'Exp. q-value: Combined', 'Coverage [%]']
            self.result['protein'] = self.raw_result['protein'][protein_cols]

            # peptide: add a new Mod_sequence column
            peptide_cols = ['Identifier', 'Mod_sequence', 'Qvality PEP',
                            'Qvality q-value', 'Master Protein Accessions',
                            '# Missed Cleavages']
            peptide = self.raw_result['peptide']
            peptide['Modifications'] = peptide['Modifications'].fillna('')
            peptide['Mod_sequence'] = peptide['Sequence'] +'.'+ peptide['Modifications']
            self.result['peptide'] = peptide[peptide_cols]

            # PSM: add a new Mod_sequence column
            PSM_cols = ['Identifier', 'Mod_sequence','Spectrum File', 'First Scan',
                        'Precursor Abundance', 'Intensity','Percolator q-Value', 'Percolator PEP',
                        'Master Scan(s)', 'Isolation Interference [%]',
                        'RT [min]','Master Protein Accessions']
            PSM = self.raw_result['PSM']
            PSM['Modifications'] = PSM['Modifications'].fillna('')
            PSM['Mod_sequence'] = PSM['Sequence'] +'.'+ PSM['Modifications']
            self.result['PSM'] = PSM[PSM_cols]
        print('Read data at {}.\nProtein table {}.\nPeptide table {}.\n'
              'PSM table {}.'.format(self.filename,
                                     self.result['protein'].shape,
                                     self.result['peptide'].shape,
                                     self.result['PSM'].shape))
        return self


def main():
    # a brief tester
    result1 = Search_Result(filename='test/two_files-PD_output/19-149_HM0522_GM-1', identifier='overnight',
                        search_workflow='ID_PD')
    result1 = result1.read_result()
    result2 = Search_Result(filename='test/two_files-PD_output/19-149_HM0522_GM-2', identifier='2x digestion',
                        search_workflow='ID_PD')
    result2 = result2.read_result()
    # print(df1.raw_result['protein'].shape)
    # print(df1.result['peptide'].shape)
    # df1.result['peptide']

    # print(df1.raw_result['PSM'].columns)

    # print(df1.result['PSM'])


if __name__ == '__main__':
    main()
