import itertools
import datetime
import argparse
import pandas as pd
import os
import re

from tol.sources.portal import portal
from tol.core import DataSourceFilter

def query_portal(plates, verbose):

    prtl = portal()
    f = DataSourceFilter()
    f.and_ = {'sts_rackid': {'in_list': {'value': plates}}}
    samples = prtl.get_list('sample', object_filters=f)
    sample_data = {}
    for i, sample in enumerate(samples):
        # print(sample.id)
        # print(sample.attributes)
        # raise Exception('debug')
        uid = sample.id
        # workaround for deleted entries
        # if sample.sts_sampleset_id is None:
        #     if verbose:
        #         print('excluding', uid, sample.sts_rackid, sample.sts_tubeid)
        #     continue
        if verbose:
            print(uid, sample.sts_rackid, sample.sts_tubeid)
        
        sample_data[uid] = {}
        sample_data[uid]['plate_id'] = sample.sts_rackid
        sample_data[uid]['well_id'] = sample.sts_tubeid
        
        sample_data[uid]['specimen_id'] = sample.sts_specimen.id
        sample_data[uid]['cohort'] = sample.sts_gal_abbreviation
        sample_data[uid]['date_of_sample_collection'] = sample.sts_col_date
        if hasattr(sample.sts_species, 'id'):
            sample_data[uid]['taxon_id'] = sample.sts_species.id
            # sample.sts_species.sts_scientific_name seems unavailable
            if sample.sts_species.id == '32644':
                sample_data[uid]['common_name'] = 'unidentified'
            elif sample.sts_species.id == '2582415':
                sample_data[uid]['common_name'] = 'blank sample'
            else:
                raise ValueError(f'taxon id {sample.sts_species.id} not expected for BIOSCAN samples')
        else:
            print(f'no species information for {sample.sts_specimen.id}')
            sample_data[uid]['taxon_id'] = '32644'
            sample_data[uid]['common_name'] = 'unidentified'

    df = pd.DataFrame(sample_data).T

    if len(df) == 0:
        df = pd.DataFrame(columns='plate_id well_id specimen_id cohort date_of_sample_collection taxon_id common_name'.split())
    
    missing_plates = set(plates) - set(df.plate_id.unique())
    extra_plates = set(df.plate_id.unique()) - set(plates)

    if len(missing_plates) > 0:
        print(f'ERROR: could not locate plates in ToL Portal {sorted(missing_plates)}')
    if len(extra_plates) > 0:
        print(f'ERROR: data extracted for extra plates {sorted(extra_plates)}')
    if len(missing_plates) == 0 and len(extra_plates) == 0:
        print('All plates were found in ToL Portal')

    return df, missing_plates

def finalise_table(df, plates, is_lysate):

    row_id = list('ABCDEFGH')
    col_id = range(1,13)
    expected_wells = [r + str(c) for (c, r) in itertools.product(col_id, row_id)]

    # sort values by plate and well
    df['plate_id'] = df['plate_id'].astype("category").cat.set_categories(plates)
    df['well_id'] = df['well_id'].astype("category").cat.set_categories(expected_wells)
    df = df.sort_values(by=['plate_id', 'well_id']).reset_index(drop=True)

    # check for plate completeness
    samples_per_plate = df['plate_id'].value_counts()
    incomplete_plates = samples_per_plate[samples_per_plate != 96]
    if len(incomplete_plates) > 0:
        print(f'ERROR: incomplete plates found: {incomplete_plates}')

    # only collection year needed. This only works expecting YYYY-MM-DD format
    # blank samples will have collection year of the previous sample
    df['date_of_sample_collection'] = df['date_of_sample_collection'].ffill().bfill().astype(str).str.split('-').str.get(0)

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

def add_sts_meta(df, plates, sts_fn):

    print(f'adding STS metadata from {sts_fn}')
    sts_df = pd.read_excel(sts_fn, sheet_name='Metadata Entry')

    added_plates_df = sts_df[sts_df.RACK_OR_PLATE_ID.isin(plates)].copy()
    if added_plates_df.shape[0] == 0:
        print(f'no plates added from {sts_fn}')
    else:
        added_plates_df['plate_id'] = added_plates_df['RACK_OR_PLATE_ID']
        added_plates_df['well_id'] = added_plates_df['TUBE_OR_WELL_ID']
        added_plates_df['specimen_id'] = added_plates_df['SPECIMEN_ID']
        cohort = sts_fn.split('/')[-1].split('.')[0].split('_')[0]
        added_plates_df['cohort'] = cohort
        added_plates_df['date_of_sample_collection'] = added_plates_df['DATE_OF_COLLECTION']
        added_plates_df['taxon_id'] = added_plates_df['TAXON_ID']
        added_plates_df['common_name'] = added_plates_df['SCIENTIFIC_NAME']

        added_plates = added_plates_df.plate_id.unique()
        print(f'added plates {", ".join(added_plates)}')

        df = pd.concat([
            df,
            added_plates_df[[
                'plate_id',
                'well_id',
                'specimen_id',
                'cohort',
                'date_of_sample_collection',
                'taxon_id',
                'common_name'
                ]]
        ]).reset_index(drop=True)

        missing_plates = [plate for plate in plates if plate not in added_plates]

    return df, missing_plates


def main():
    
    parser = argparse.ArgumentParser("Generate bioscan sciops manifest given a list of plate IDs. "
        "By default, specimen (LILYS) manifest is generated.")
    parser.add_argument('-p', '--plates', help='File listing plates, one per line', required=True)
    parser.add_argument('-s', '--sts_manifests', 
        help='STS manifest used to add sample info prior to querying portal', 
        action='append')
    parser.add_argument(
        '-o', '--outfile', 
        help='Output file. Default: results/[lilys/lbsn]_YYYYMMDD.xlsx', 
        default=None)
    parser.add_argument('-l', '--lysate', help='Generate manifest for lysate (LBSN) plates instead: '
                        'add positive control at G12', action='store_true')
    parser.add_argument('-v', '--verbose', help='Include sample-level messages', 
                        action='store_true', default=False)

    args = parser.parse_args()

    if args.outfile is None:
        today = datetime.datetime.now().strftime('%Y%m%d')
        if args.lysate:
            args.outfile = f'results/lbsn_{today}.xlsx'
        else:
            args.outfile = f'results/lilys_{today}.xlsx'
    assert args.outfile.endswith('tsv') or args.outfile.endswith('xlsx'), 'can only write to ".tsv" or ".xlsx" file'

    plates = []
    with open(args.plates) as f:
        for line in f:
            plate = line.strip()
            assert ' ' not in plate, f'plate "{plate}" contains space character, aborting'
            if len(plate) > 0:
                plates.append(plate)

    df = pd.DataFrame()

    if args.sts_manifests is not None:
        for sts_fn in args.sts_manifests:
            assert os.path.isfile(sts_fn)
            cohort = sts_fn.split('/')[-1].split('.')[0].split('_')[0]
            assert re.match(r'^[A-Z]{4}$', cohort), (
                f'we need STS manifest filename to start with four-letter partner code '
                f'followed by underscore, found {cohort} instead'
                )
            df, plates = add_sts_meta(df, plates, sts_fn)
            # sts_sampleset = sts_fn.split('/')[-1].split('.')[0]
            # assert re.match(r'^[A-Z]{4}_[0-9]{6}$', sts_sampleset), (
            #     f'we need STS manifest filename to be STS sampleset ID like ABCD_123456, '
            #     f'found {sts_sampleset} instead'
            #     )

    max_plates = 104
    if len(plates) > max_plates:
        raise ValueError(
            f'Can only query ToL portal for 10,000 samples max - this is {max_plates} plates'
            )
    if len(plates) == 0:
        print('Skipping portal query as all plates added from STS manifest')
    else:
        print(f'Querying ToL Portal for {len(plates)} plates')
        portal_df, missing_plates = query_portal(plates, args.verbose)
        df = pd.concat([
            df,
            portal_df
            ]).reset_index(drop=True)

        if len(missing_plates) > 0:
            print(f'Could not find plates on portal: {sorted(missing_plates)}')
        else:
            print(f'Found all {len(plates)} plates')

    plate_type = 'lysate (LBSN)' if args.lysate else 'specimen (LILYS)'
    print(f'Adjusting tables for {plate_type} plate SciOps submission ')
    df = finalise_table(df, plates, args.lysate)
    print(f'Writing to {args.outfile}')
    if args.outfile.endswith('tsv'):
        df.to_csv(args.outfile, sep='\t', index=False)
    elif args.outfile.endswith('xlsx'):
        df.to_excel(args.outfile, index=False)


if __name__ == '__main__':
    main()