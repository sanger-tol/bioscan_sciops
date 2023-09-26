import argparse
import psycopg2
import pandas as pd
from configparser import ConfigParser


def read_config(filename, section='postgresql'):

    # create a parser
    parser = ConfigParser()
    # read config file
    parser.read(filename)

    # get section, default to postgresql
    db = {}
    if parser.has_section(section):
        params = parser.items(section)
        for param in params:
            db[param[0]] = param[1]
    else:
        raise Exception('Section {0} not found in the {1} file'.format(section, filename))

    return db


def query_sts(plates, config):

    conn = None
    try:
        # read connection parameters
        params = read_config(config)

        # connect to the PostgreSQL server
        print('Connecting to the PostgreSQL database...')
        conn = psycopg2.connect(**params)
        
        # create a cursor
        cur = conn.cursor()
        
        # execute a statement
        print('PostgreSQL database version:')
        cur.execute('SELECT version()')

        # display the PostgreSQL database server version
        db_version = cur.fetchone()
        print(db_version)

        # template query to check table columns
        # dev_query = ('''
        #     SELECT
        #         attname AS colname,
        #         pg_catalog.format_type(atttypid, atttypmod) AS coltype
        #     FROM
        #         pg_catalog.pg_attribute
        #     WHERE
        #         attnum > 0 AND
        #         NOT attisdropped AND
        #         attrelid = 'gal'::regclass
        #     ORDER BY
        #         attnum ASC
        #     ''')

        # Execute a PostgreSQL query
        
        # collection country not properly recorded - sample.collection_country_id is always None
        # sample.loc_id = location.location_id -> location.location - does not include country name

        query = (f'''
            select sample.manifest_id, sample.series, sample.tubeid, sample.specimenid, gal.abbreviation, sample.col_date, 
                species.taxonid, species.scientific_name, sample.rackid
            from sample, sample_species, species, gal
            where sample.sample_id = sample_species.sample_id
                and sample_species.species_id = species.species_id
                and sample.gal_id = gal.gal_id
                and rackid in {str(plates).replace('[','(').replace(']',')')}
            ''')

        cur.execute(query)

        # Fetch and print the results
        rows = cur.fetchall()
        # cur.execute(dev_query)
        # dev_rows = cur.fetchall()
        # for row in dev_rows:
        #     print(row)
        colnames = ['manifest', 'series', 'well_id', 'specimen_id', 'cohort', 'date_of_sample_collection', 'taxon_id', 'common_name', 'sample_description']
        df = pd.DataFrame(rows, columns=colnames)
        df['series'] = df['series'].astype(int)
        df = df.sort_values(by=['manifest','series']).reset_index(drop=True)
        print(df)
       
        # close the communication with the PostgreSQL
        cur.close()
    except (Exception, psycopg2.DatabaseError) as error:
        print(error)
    finally:
        if conn is not None:
            conn.close()
            print('Database connection closed.')

    return df

def finalise_table(df):

    # auto-fill
    df['country_of_origin'] = 'United Kingdom'
    df['retention_instruction'] = 'Return to customer after 2 years'

    # mark up controls
    df['bioscan_supplier_sample_name'] = df['specimen_id']
    df['bioscan_control_type'] = ''
    pos_controls = ((df.common_name == 'blank sample') & (df.well_id == 'G12'))
    df.loc[pos_controls, 'bioscan_supplier_sample_name'] = 'CONTROL_POS_' + df.loc[pos_controls, 'specimen_id']
    df.loc[pos_controls, 'bioscan_control_type'] = 'pcr positive'
    neg_controls = ((df.common_name == 'blank sample') & (df.well_id != 'G12'))
    df.loc[neg_controls, 'bioscan_supplier_sample_name'] = 'CONTROL_NEG_LYSATE_' + df.loc[neg_controls, 'specimen_id']
    df.loc[neg_controls, 'bioscan_control_type'] = 'lysate negative'

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
    
    parser = argparse.ArgumentParser("Generate bioscan sciops manifest given a list of plate IDs")
    parser.add_argument('-p', '--plates', help='File listing plates, one per line', required=True)
    parser.add_argument('-c', '--config', help='Config file with STS database credentials', default='../sts_config.ini')
    # parser.add_argument('-s', '--stats', help='DADA2 stats tsv file', required=True)
    parser.add_argument('-o', '--outfile', help='Output file. Default: out.tsv', default='out.tsv')
    parser.add_argument('-v', '--verbose', 
                        help='Include INFO level log messages', action='store_true')

    args = parser.parse_args()
    
    assert args.outfile.endswith('tsv'), 'can only write to ".tsv" file'
    
    plates = []
    with open(args.plates) as f:
        for line in f:
            plate = line.strip()
            if len(plate) > 0:
                plates.append(plate)

    print(f'Querying STS for {len(plates)} plates')
    df = query_sts(plates, args.config, args.verbose)
    df = finalise_table(df)
    print(f'Writing to {args.outfile}')
    df.to_csv(args.outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()