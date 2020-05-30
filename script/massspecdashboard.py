import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# from matplotlib_venn import venn3

class Search_Result():
    def __init__(self, path='', identifier='',
                 search_workflow=''):
        """use resultname.path, resultname.identifier, etc.

        could extend more using a method for i/o.
        also could save in class.
        """
        self.path = path  # implies where to read
        self.identifier = identifier
        self.search_workflow = search_workflow  # ID/Quan_search engine (top-down hierarchy)
        if search_workflow == 'ID_PD':
            self.read_result_ID_PD()
        elif search_workflow == 'ID_MQ':
            self.read_result_ID_MQ()

    def read_result_ID_PD(self):
        # output: self.raw_result, self.result
        # a dict {protein: protein_result, peptide: peptide_result, PSM: PSM_result}
        self.raw_result = {}
        if self.search_workflow == 'ID_PD':
            suffix = {'PSM': self.path + '_PSMs.txt',
                      'peptide': self.path + '_PeptideGroups.txt',
                      'protein': self.path + '_Proteins.txt'}
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
            PSM_cols = [
                'Identifier', 'Mod_sequence','Spectrum File', 'First Scan',
                'Precursor Abundance', 'Intensity','Percolator q-Value',
                'Percolator PEP', 'Master Scan(s)', 'Isolation Interference [%]',
                'RT [min]','Master Protein Accessions'
            ]
            PSM = self.raw_result['PSM']
            PSM['Modifications'] = PSM['Modifications'].fillna('')
            PSM['Mod_sequence'] = PSM['Sequence'] +'.'+ PSM['Modifications']
            self.result['PSM'] = PSM[PSM_cols]

        print('Read data ID_PD of identifier {} at {}.\nProtein table {}.\nPeptide table {}.\n'
              'PSM table {}.'.format(self.identifier,
                                     self.path,
                                     self.result['protein'].shape,
                                     self.result['peptide'].shape,
                                     self.result['PSM'].shape))
        self.protein = self.result['protein']
        self.peptide = self.result['peptide']
        self.PSM = self.result['PSM']
        return self

    def read_result_ID_MQ(self):
        # output: self.raw_result, self.result
        # a dict {protein: protein_result, peptide: peptide_result, PSM: PSM_result}
        # TODO(mofhu): rename columns to same format as PD
        self.raw_result = {}
        if self.search_workflow == 'ID_MQ':
            suffix = {'PSM': self.path + '/evidence.txt',
                      'peptide': self.path + '/peptides.txt',
                      'protein': self.path + '/proteinGroups.txt',
                      'summary': self.path + '/summary.txt'}
            for key in suffix:
                self.raw_result[key] = pd.read_csv(suffix[key], sep='\t', dtype={'Experiment':str})
                self.raw_result[key]['Identifier'] = self.identifier

        # raw result to final result
        self.result = {}
        if self.search_workflow == 'ID_MQ':
            # select columns
            # protein: no change, only filter columns
            protein = self.raw_result['protein']
            protein = protein.rename(columns={
                'Majority protein IDs': 'Accession'
            })
            protein_cols = ['Identifier', 'Accession',
                            'Sequence coverage [%]',
                            'Razor + unique peptides']
            self.result['protein'] = protein[protein_cols]
            # peptide: add a new Mod_sequence column
            peptide = self.raw_result['peptide']
            peptide = peptide.rename(columns={
                'Sequence': 'Mod_sequence'
            })
            peptide_cols = ['Identifier', 'Mod_sequence', 'Missed cleavages',
                            'PEP', 'MS/MS Count']
            self.result['peptide'] = peptide[peptide_cols]
            # PSM: add a new Mod_sequence column
            PSM = self.raw_result['PSM']
            PSM = PSM.rename(columns={
                'Sequence': 'Mod_sequence'
            })
            PSM_cols = [
                'Identifier', 'Mod_sequence',
                'PIF', 'Fraction of total spectrum',
                'Base peak fraction', 'PEP',
            ]
            self.result['PSM'] = PSM[PSM_cols]
        # protein number is not filtered, as we concentrate on peptides and PSM first, this is fine.
        print('Read data ID_MQ of identifier {} at {}.\n'
              'Peptide table {}.\n'
              'PSM table {}.\nProtein table {}.'.format(self.identifier,
                                     self.path,
                                     self.result['protein'].shape[0],
                                     self.result['peptide'].shape[0],
                                     self.result['PSM'].shape[0]))
        self.protein = self.result['protein']
        self.peptide = self.result['peptide']
        self.PSM = self.result['PSM']
        return self

def filter_protein_count(protein_db, threshold=None, search_workflow=''):
    # protein filter
    if search_workflow == 'ID_PD':
        protein_fdr = protein_db[protein_db['Exp. q-value: Combined'] < 0.01]
        protein_master = protein_fdr.loc[protein_fdr['Master'].str.match('IsMasterProtein$')]  # use regex to filter IsMasterProteinCandidate
        pivot = pd.pivot_table(protein_master, values='Accession',
                            index='Identifier', aggfunc='count')
    elif search_workflow == 'ID_MQ':
        pivot = pd.pivot_table(protein_db, values='Accession',
                            index='Identifier', aggfunc='count')
    # flatten pivot table to normal df
    flatten = pd.DataFrame(pivot.to_records())
    if threshold is not None:
        flatten = flatten.append(
            {
                'Identifier': "THRESHOLD",
                'Accession': threshold,
            },
            ignore_index=True)
    print('protein count after FDR cutoff:')
    print(flatten)
    return flatten


def filter_peptide_count(peptide_db, threshold=None, search_workflow=''):
    # peptide filter
    pivot = pd.pivot_table(peptide_db, values='Mod_sequence',
                           index='Identifier', aggfunc='count')
    # flatten pivot table to normal df
    flatten = pd.DataFrame(pivot.to_records())
    if threshold is not None:
        flatten = flatten.append(
            {
                'Identifier': "THRESHOLD",
                'Mod_sequence': threshold,
            },
            ignore_index=True)
    print('peptide count:')
    print(flatten)
    return flatten


def filter_PSM_count(PSM_db, threshold=None, search_workflow=''):
    # peptide filter
    pivot = pd.pivot_table(PSM_db, values='Mod_sequence',
                           index='Identifier', aggfunc='count')
    # flatten pivot table to normal df
    flatten = pd.DataFrame(pivot.to_records())
    if threshold is not None:
        flatten = flatten.append(
            {
                'Identifier': "THRESHOLD",
                'Mod_sequence': threshold,
            },
            ignore_index=True)
    print('PSM count:')
    print(flatten)
    return flatten


def barplot(data, x, y, path='', output='png'):
    # set seaborn and annotation codes
    plt.figure(figsize=(6, 8))
    sns.set(style='whitegrid')
    splot = sns.barplot(data=data, x=x, y=y)
    for p in splot.patches:
        splot.annotate(format(p.get_height(), ), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'center', xytext = (0, 10), textcoords = 'offset points')
    if output == 'png':
        plt.title(path)
        plt.savefig(path+'.png')
    else:
        output.savefig()


def main():
    """Main function for Mass Spec DashBoard.
    input: config file format YAML:
        - identifier: file name (no suffix)
    output:
        - a TXT file of meta data and result DataFrame
        - several PNG figures of protein, peptide, and PSM
        - a PDF file of three level figures
    """
    import yaml
    parser = argparse.ArgumentParser()
    # one or more path
    parser.add_argument('config_file', help='config file in yaml format.')
    args = parser.parse_args()
    config_file = args.config_file
    config = yaml.load(open(config_file).read(), Loader=yaml.CLoader)

    from time import localtime, asctime
    file_out = open('massspecdashboard.log', 'w')
    file_out.write('MassSpecDashBoard' + '\n')
    file_out.write('Timestamp: ' + asctime(localtime()) + '\n')
    file_out.write('Config file:\n')
    file_out.write(open(config_file).read())

    result = {'protein':[],
              'peptide':[],
              'PSM':[]}
    workflow = config['WORKFLOW']
    if 'THRESHOLD' in config:
        threshold = config['THRESHOLD']
    else:
        # no threshold: init with None
        threshold = {
            'protein': None,
            'peptide': None,
            'PSM': None,
        }
    for identifier in config:
        if identifier == 'THRESHOLD' or identifier == 'WORKFLOW':
            continue
        path = config[identifier]
        tables = Search_Result(path=path, identifier=identifier,
                                search_workflow=workflow)
        result['protein'].append(tables.protein)
        result['peptide'].append(tables.peptide)
        result['PSM'].append(tables.PSM)
    # protein result
    # TODO(mofhu): add PD protein filter in class, not in function?
    protein_result = filter_protein_count(pd.concat(result['protein']), threshold=threshold['protein'], search_workflow=workflow)
    file_out.write('---\nResult:\n')
    file_out.write('Protein:\n' + protein_result.to_string(index=False) + '\n')
    barplot(data=protein_result, x='Identifier', y='Accession', path='protein')
    # peptide result
    peptide_result = filter_peptide_count(pd.concat(result['peptide']), threshold=threshold['peptide'], search_workflow=workflow)
    file_out.write('Peptide:\n' + peptide_result.to_string(index=False) + '\n')
    barplot(data=peptide_result, x='Identifier', y='Mod_sequence', path='peptide')
    # PSM result
    PSM_result = filter_PSM_count(pd.concat(result['PSM']), threshold=threshold['PSM'], search_workflow = workflow)
    file_out.write('PSM:\n' + PSM_result.to_string(index=False) + '\n')
    barplot(data=PSM_result, x='Identifier', y='Mod_sequence', path='PSM')
    # pdf single report
    from matplotlib.backends.backend_pdf import PdfPages
    # TODO(mofhu): update pdf export to a better output clarity
    pp = PdfPages('massspecdashboard-report.pdf')

    # TODO(mofhu): add a plot name
    barplot(data=protein_result, x='Identifier', y='Accession', path='protein', output=pp)
    barplot(data=peptide_result, x='Identifier', y='Mod_sequence', path='peptide', output=pp)
    barplot(data=PSM_result, x='Identifier', y='Mod_sequence', path='PSM', output=pp)
    pp.close()

if __name__ == '__main__':
    main()
