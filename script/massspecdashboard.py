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
        self.read_result()

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
            protein_cols = ['Identifier', 'Accession', 'Master',
                            'Exp. q-value: Combined', 'Coverage [%]']
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
        self.protein = self.result['protein']
        self.peptide = self.result['peptide']
        self.PSM = self.result['PSM']
        return self


def filter_protein_count(protein_db):
    # protein filter
    protein_fdr = protein_db[protein_db['Exp. q-value: Combined'] < 0.01]
    protein_master = protein_fdr.loc[protein_fdr['Master'].str.match('IsMasterProtein')]
    pivot = pd.pivot_table(protein_master, values='Accession',
                           index='Identifier', aggfunc='count')
    # flatten pivot table to normal df
    flatten = pd.DataFrame(pivot.to_records())
    print('protein count after FDR cutoff:')
    print(flatten)
    return flatten


def filter_peptide_count(peptide_db):
    # peptide filter
    pivot = pd.pivot_table(peptide_db, values='Mod_sequence',
                           index='Identifier', aggfunc='count')
    # flatten pivot table to normal df
    flatten = pd.DataFrame(pivot.to_records())
    print('peptide count:')
    print(flatten)
    return flatten


def filter_PSM_count(PSM_db):
    # peptide filter
    pivot = pd.pivot_table(PSM_db, values='Mod_sequence',
                           index='Identifier', aggfunc='count')
    # flatten pivot table to normal df
    flatten = pd.DataFrame(pivot.to_records())
    print('PSM count:')
    print(flatten)
    return flatten


def barplot(data, x, y, filename=''):
    # set seaborn and annotation codes
    plt.figure(figsize=(6, 8))
    sns.set(style='whitegrid')
    splot = sns.barplot(data=data, x=x, y=y)
    for p in splot.patches:
        splot.annotate(format(p.get_height(), ), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'center', xytext = (0, 10), textcoords = 'offset points')
    plt.title(filename)
    plt.savefig(filename+'.png')


def main():
    # a brief tester
    result1 = Search_Result(filename='test/two_files-PD_output/19-149_HM0522_GM-1', identifier='1x digestion',
                        search_workflow='ID_PD')
    # result1 = result1.read_result()
    result2 = Search_Result(filename='test/two_files-PD_output/19-149_HM0522_GM-2', identifier='2x digestion',
                        search_workflow='ID_PD')
    # result2 = result2.read_result()
    result3 = Search_Result(filename='test/two_files-PD_output/19-149_HM0522_GM-2', identifier='2x digestion-rep',
                        search_workflow='ID_PD')
    # result3 = result3.read_result()
    protein_result = filter_protein_count(pd.concat([result1.protein, result2.protein,result3.protein]))
    barplot(data=protein_result, x='Identifier', y='Accession', filename='protein')
    # peptide_result
    peptide_result = filter_peptide_count(pd.concat([result1.peptide, result2.peptide, result3.peptide]))
    barplot(data=peptide_result, x='Identifier', y='Mod_sequence', filename='peptide')
    # PSM_result
    PSM_result = filter_PSM_count(pd.concat([result1.PSM, result2.PSM, result3.PSM]))
    barplot(data=PSM_result, x='Identifier', y='Mod_sequence', filename='PSM')

if __name__ == '__main__':
    main()
