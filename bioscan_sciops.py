import itertools
import argparse
import pandas as pd

from tol.sources.portal import portal
from tol.core import DataSourceFilter


def query_portal(plates, verbose):

    prtl = portal()
    f = DataSourceFilter()
    f.in_list = {'sts_rackid': plates}
    samples = prtl.get_list('sample', object_filters=f)
    sample_data = {}
    for i, sample in enumerate(samples):
        uid = sample.uid
        # workaround for deleted entries
        if sample.sts_sampleset_id is None:
            if verbose:
                print('excluding', uid, sample.sts_rackid, sample.sts_tubeid)
            continue
        if verbose:
            print(uid, sample.sts_rackid, sample.sts_tubeid)
        
        sample_data[uid] = {}
        sample_data[uid]['plate_id'] = sample.sts_rackid
        sample_data[uid]['well_id'] = sample.sts_tubeid
        
        sample_data[uid]['specimen_id'] = sample.sts_specimen.id
        sample_data[uid]['cohort'] = sample.sts_gal_abbreviation
        sample_data[uid]['date_of_sample_collection'] = sample.sts_col_date
        sample_data[uid]['taxon_id'] = sample.sts_species.id
        # sample.sts_species.sts_scientific_name seems unavailable
        if sample.sts_species.id == '32644':
            sample_data[uid]['common_name'] = 'unidentified'
        elif sample.sts_species.id == '2582415':
            sample_data[uid]['common_name'] = 'blank sample'
        else:
            raise ValueError(f'taxon id {sample.sts_species.id} not expected for BIOSCAN samples')

    df = pd.DataFrame(sample_data).T
    
    missing_plates = set(plates) - set(df.plate_id.unique())
    extra_plates = set(df.plate_id.unique()) - set(plates)

    if len(missing_plates) > 0:
        print(f'ERROR: could not locate plates in ToL Portal {missing_plates}')
    if len(extra_plates) > 0:
        print(f'ERROR: data extracted for extra plates {extra_plates}')
    if len(missing_plates) == 0 and len(extra_plates) == 0:
        print('All plates were found in ToL Portal')

    return df


def finalise_table(df, plates, is_lysate):

    row_id = list('ABCDEFGH')
    col_id = range(1,13)
    expected_wells = [r + str(c) for (c, r) in itertools.product(col_id, row_id)]

    # sort values by plate and well
    df['plate_id'] = df['plate_id'].astype("category").cat.set_categories(plates)
    df['well_id'] = df['well_id'].astype("category").cat.set_categories(expected_wells)
    df = df.sort_values(by=['plate_id', 'well_id']).reset_index(drop=True)

    # only collection year needed. This only works expecting YYYY-MM-DD format
    # blank samples will have collection year of the previous sample
    df['date_of_sample_collection'] = df['date_of_sample_collection'].ffill().astype(str).str.split('-').str.get(0)
    if df['date_of_sample_collection'].isna().any():
        print('Warning: empty time entries (NaT) were generated - please fix manually')

    # auto-fill
    df['country_of_origin'] = 'United Kingdom'
    df['retention_instruction'] = 'Return to customer after 2 years'
    df['sample_description'] = df['plate_id']

    # mark up controls
    df['bioscan_supplier_sample_name'] = df['specimen_id']
    df['bioscan_control_type'] = ''
    # pos control does not have to be (df.common_name == 'blank sample')
    # we do not modify positive controls for specimen plates
    # we do add positive control for lysis plates
    if is_lysate:
        pos_controls = (df.well_id == 'G12')
        df.loc[pos_controls, 'bioscan_supplier_sample_name'] = 'CONTROL_POS_' + df.loc[pos_controls, 'specimen_id']
        df.loc[pos_controls, 'taxon_id'] = '32644'
        df.loc[pos_controls, 'common_name'] = 'unidentified'
        # we do not touch sciops lims control type in case of lysate
        # df.loc[pos_controls, 'bioscan_control_type'] = 'pcr positive'
    neg_controls = ((df.common_name == 'blank sample') & (df.well_id != 'G12'))
    df.loc[neg_controls, 'bioscan_supplier_sample_name'] = 'CONTROL_NEG_LYSATE_' + df.loc[neg_controls, 'specimen_id']
    # only mark H12 as lysate negative in control type
    df.loc[neg_controls & (df.well_id == 'H12'), 'bioscan_control_type'] = 'lysate negative'

    # sanity check taxonomy
    expected_taxa = ['unidentified','blank sample']
    if not df.common_name.isin(expected_taxa).all():
        unexpected_taxa = set(df.common_name.unique) - set(expected_taxa)
        unexpected_taxa_samples = df.common_name.isin(unexpected_taxa).bioscan_supplier_sample_name.to_list()
        print('ERROR: found unexpected taxa {unexpected_taxa} for samples {unexpected_taxa_samples}')

    # reorder columns
    out_df = df[[
        'bioscan_supplier_sample_name',
        'retention_instruction',
        'cohort',
        'country_of_origin',
        'date_of_sample_collection',
        'taxon_id',
        'common_name',
        'sample_description',
        'bioscan_control_type'
    ]].copy()

    out_df.columns = out_df.columns.str.replace('_',' ').str.upper()    

    return out_df


def main():
    
    parser = argparse.ArgumentParser("Generate bioscan sciops manifest given a list of plate IDs. "
        "By default, specimen (LILYS) manifest is generated.")
    parser.add_argument('-p', '--plates', help='File listing plates, one per line', required=True)
    parser.add_argument('-o', '--outfile', help='Output file. Default: out.tsv', default='out.tsv')
    parser.add_argument('-l', '--lysate', help='Generate manifest for lysate (LBSN) plates instead: '
                        'add positive control at G12', action='store_true')
    parser.add_argument('-v', '--verbose', help='Include sample-level messages', 
                        action='store_true', default=False)

    args = parser.parse_args()
    
    assert args.outfile.endswith('tsv'), 'can only write to ".tsv" file'
    
    plates = []
    with open(args.plates) as f:
        for line in f:
            plate = line.strip()
            assert ' ' not in plate, f'plate "{plate}" contains space character, aborting'
            if len(plate) > 0:
                plates.append(plate)

    max_plates = 104
    if len(plates) > max_plates:
        raise ValueError(
            f'Can only query ToL portal for 10,000 samples max - this is {max_plates} plates'
            )

    print(f'Querying ToL Portal for {len(plates)} plates')
    df = query_portal(plates, args.verbose)
    plate_type = 'lysate (LBSN)' if args.lysate else 'specimen (LILYS)'
    print(f'Adjusting tables for {plate_type} plate SciOps submission ')
    df = finalise_table(df, plates, args.lysate)
    print(f'Writing to {args.outfile}')
    df.to_csv(args.outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()