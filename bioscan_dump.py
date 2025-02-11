import datetime
import argparse
import pandas as pd
import os

from tol.sources.portal import portal
from tol.core import DataSourceFilter

def query_portal(samplesets):

    prtl = portal()
    f = DataSourceFilter()
    os.makedirs('results/sampleset_dumps', exist_ok=True)
    for sampleset in samplesets:
        sampleset_fn = f'results/sampleset_dumps/{sampleset}.dump.csv'
        if os.path.isfile(sampleset_fn):
            print(f'dump exists at {sampleset_fn}, skipping')
            continue
        print(f'fetching {sampleset} from portal')
        f.and_ = {'sts_project': {'eq': {'value': 'BIOSCAN'}}, 'sts_sampleset.id': {'eq': {'value': sampleset}}}
        samples = prtl.get_list('sample', object_filters=f)
        print('unpacking...')
        sampletadata = [sample._CoreDataObject__attributes for sample in samples]
        if len(sampletadata) == 0:
            raise ValueError(f'no samples found for {sampleset} in bioscan')
        pd.DataFrame(sampletadata).to_csv(sampleset_fn, index=False)
        print(f'saved to {sampleset_fn}')
        

def combine_dumps(samplesets):

    sampleset_dfs = []
    for sampleset in samplesets:
        sampleset_fn = f'results/sampleset_dumps/{sampleset}.dump.csv'
        sampleset_df = pd.read_csv(sampleset_fn)
        sampleset_dfs.append(sampleset_df)

    return pd.concat(sampleset_dfs)

def main():
    
    parser = argparse.ArgumentParser("Generate bioscan data dump given a list of manifest IDs. ")
    parser.add_argument('samplesets', help='File listing samplesets, one per line')

    args = parser.parse_args()

    samplesets = []
    with open(args.samplesets) as f:
        for line in f:
            sampleset = line.strip()
            assert ' ' not in sampleset, f'plate "{sampleset}" contains space character, aborting'
            if len(sampleset) > 0: # skipping empty lines
                samplesets.append(sampleset)
    print(f'found {len(sampleset)} samplesets, proceeding to portal query')

    query_portal(samplesets)

    print('combining datasets currently disabled, need to subset columns first')
    # df = combine_dumps(samplesets)

    # now = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    # outfile = f'results/bioscan_dump_{now}.csv'
    # df.to_csv(outfile, index=False)

    print('all done')


if __name__ == '__main__':
    main()