from __future__ import print_function
# from __future__ import unicode_literals


def entrez_db_list(email):

    '''
    Returns:
        A list of Entrez databases sorted alphabetically.
    '''

    from Bio import Entrez
    Entrez.email = email
    handle = Entrez.einfo()
    record = Entrez.read(handle)
    dbs = record['DbList']
    dbs.sort()
    return dbs


def esearch(esearch_terms, db, email):

    '''
    Perform Entrez ESearch by term.

    Args:
        esearch_terms: One or more search terms that use esearch syntax
            http://www.ncbi.nlm.nih.gov/books/NBK3837/
                #EntrezHelp.Entrez_Searching_Options

        db: Entrez database name. Use entrez_db_list to get a current list of
            available databases.
            pubmed, protein, nuccore, nucleotide, nucgss, nucest, structure,
            genome, assembly, gcassembly, genomeprj, bioproject, biosample,
            biosystems, blastdbinfo, books, cdd, clone, gap, gapplus, dbvar,
            epigenomics, gene, gds, geo, geoprofiles, homologene, journals,
            mesh, ncbisearch, nlmcatalog, omia, omim, pmc, popset, probe,
            proteinclusters, pcassay, pccompound, pcsubstance, pubmedhealth,
            seqannot, snp, sra, taxonomy, toolkit, toolkitall, unigene, unists,
            gencoll

        email: An email address at which the user can be reached. To make use
            of NCBI's E-utilities, NCBI strongly recommends you to specify your
            email address with each request. In case of excessive usage of the
            E-utilities, NCBI will attempt to contact a user at the email
            address provided before blocking access to the E-utilities.

    Returns:
        A set of UID's.
    '''

    import time
    from Bio import Entrez
    Entrez.email = email
    if isinstance(esearch_terms, basestring):
        esearch_terms = [esearch_terms]
    retmax = None
    uid_set = []
    i = 0
    while True:
        try:
            for term in esearch_terms:
                handle = Entrez.egquery(term=term)
                record = Entrez.read(handle)
                for row in record['eGQueryResult']:
                    if row['DbName'] == db:
                        retmax = int(row['Count'])
                        break
                handle = Entrez.esearch(db=db, term=term, retmax=retmax)
                record = Entrez.read(handle)
                uid_set = uid_set + record['IdList']
            uid_set = set(uid_set)
        except:
            print('    HTTP problem, retrying...')
            i = i + 1
            if i == 10:
                break
            time.sleep(2 * i)
            continue
        break

    return uid_set


def download_ncbi_records(uids, db, entrez_email, rettype, large_batch_size=2000, small_batch_size=500, genbank_mode=False, verbose=False):

    '''
    Will download sequence records for uids and database (db) given from NCBI.
    '''

    # import time

    from Bio import Entrez
    # from Bio import SeqIO
    import krbioio

    if isinstance(uids, set):
        uids = list(uids)

    if isinstance(uids, basestring):
        uids = [uids]

    Entrez.email = entrez_email
    # out_handle = open(file_path, 'w')
    uid_count = len(uids)

    fixed_uids = list()
    for uid in uids:
        fixed_uids.append(str(uid))
    uids = fixed_uids

    # Not sure if these should be input as function arguments.
    # large_batch_size = 2000
    # small_batch_size = 500

    # Perhaps these may be function arguments?
    # rettype = 'gb'
    retmode = 'xml'

    if genbank_mode:
        retmode = 'text'

    # missing_uids = set()
    # small_batch_forced = None

    all_results = list()

    for uid_start in range(0, uid_count, large_batch_size):
        # if uid_start < 1:
        #     continue
        while True:
            # ##
            # downloaded_uids = set()
            # to_download_uids = set()
            # ##
            temp_records = list()
            n_rec_to_download = 0
            rec_downloaded = 0
            try:
                uid_end = min(uid_count, uid_start + large_batch_size)
                if verbose:
                    print('Downloading records %i to %i of %i.'
                          % (uid_start + 1, uid_end, uid_count))
                small_batch = uids[uid_start:uid_end]
                # if small_batch_forced is not None:
                #     small_batch = small_batch_forced
                if len(small_batch) == 0:
                    break
                ##
                # to_download_uids |= set(small_batch)
                ##
                small_batch_count = len(small_batch)
                small_batch_text = ','.join(small_batch)
                epost = Entrez.read(Entrez.epost(db, id=small_batch_text))
                webenv = epost['WebEnv']
                query_key = epost['QueryKey']

                temp_records = list()

                for start in range(0, small_batch_count, small_batch_size):
                    end = min(small_batch_count, start + small_batch_size)
                    if verbose:
                        print ('  Going to download record %i to %i of %i.'
                               % (start + 1, end, small_batch_count))

                    # for i, j in enumerate(range(start, end)):
                    #     print(i, small_batch[j])

                    fetch_handle = Entrez.efetch(
                        db=db, rettype=rettype, retmode=retmode,
                        retstart=start, retmax=small_batch_size,
                        webenv=webenv, query_key=query_key)

                    batch_data = None
                    if retmode == 'xml':
                        batch_data = Entrez.read(fetch_handle)
                    elif genbank_mode:
                        batch_data = krbioio.read_sequence_data(fetch_handle, rettype)
                    fetch_handle.close()
                    if batch_data:
                        temp_records = temp_records + batch_data

                # n_rec_to_download = uid_end - uid_start
                n_rec_to_download = len(small_batch)
                rec_downloaded = len(temp_records)
                ##
                # for tr in temp_records:
                #     downloaded_uids.add(tr.annotations['gi'])
                ##
            except Exception as err:
                # print(rec_downloaded, n_rec_to_download)
                # print(len(downloaded_uids), len(to_download_uids))
                # print(downloaded_uids - to_download_uids)
                # print(to_download_uids - downloaded_uids)

                print(err)

                # print('    HTTP problem, retrying...')
                # time.sleep(5)
                # continue

            # print(rec_downloaded, n_rec_to_download)
            # print(len(downloaded_uids), len(to_download_uids))
            # print(downloaded_uids - to_download_uids)
            # print(to_download_uids - downloaded_uids)

            # if rec_downloaded == n_rec_to_download:
            if verbose:
                print('    Downloaded', rec_downloaded, 'of',
                      n_rec_to_download, 'records.')
            # SeqIO.write(temp_records, out_handle, 'gb')
            all_results = all_results + temp_records
            # small_batch_forced = None
            break
            # else:
            #     missing_uids_now = to_download_uids - downloaded_uids
            #     if len(missing_uids_now & missing_uids) > 0:
            #         print('    Excluding uids:', ', '.join(list(missing_uids_now)))
            #         small_batch_forced = list(set(small_batch) - missing_uids_now)
            #         continue
            #     # print('    Downloaded', rec_downloaded, 'of',
            #     #       n_rec_to_download, 'records.')
            #     print('    Download corrupted, retrying...')
            #     print('    Missing uids:', ', '.join(list(missing_uids_now)))
            #     missing_uids |= missing_uids_now
            #     continue

    # out_handle.close()

    return all_results


def download_sequence_records(file_path, uids, db, entrez_email, verbose=True):

    '''
    Will download sequence records for uids and database (db) given from NCBI.
    '''

    from Bio import SeqIO

    records = download_ncbi_records(
        uids, db, entrez_email, rettype='gb',
        large_batch_size=2000, small_batch_size=500, genbank_mode=True, verbose=verbose)

    with open(file_path, 'w') as out_handle:
        SeqIO.write(records, out_handle, 'gb')


# def download_sequence_records(file_path, uids, db, entrez_email):

#     '''
#     Will download sequence records for uids and database (db) given from NCBI.
#     '''

#     # import time

#     from Bio import Entrez
#     from Bio import SeqIO
#     import krbioio

#     if isinstance(uids, set):
#         uids = list(uids)

#     if isinstance(uids, basestring):
#         uids = [uids]

#     Entrez.email = entrez_email
#     out_handle = open(file_path, 'w')
#     uid_count = len(uids)

#     fixed_uids = list()
#     for uid in uids:
#         fixed_uids.append(str(uid))
#     uids = fixed_uids

#     # Not sure if these should be input as function arguments.
#     large_batch_size = 2000
#     small_batch_size = 500

#     # Perhaps these may be function arguments?
#     rettype = 'gb'
#     retmode = 'text'

#     missing_uids = set()
#     small_batch_forced = None

#     for uid_start in range(0, uid_count, large_batch_size):
#         # if uid_start < 1:
#         #     continue
#         while True:
#             # ##
#             downloaded_uids = set()
#             to_download_uids = set()
#             # ##
#             try:
#                 uid_end = min(uid_count, uid_start + large_batch_size)
#                 print('Downloading records %i to %i of %i.'
#                       % (uid_start + 1, uid_end, uid_count))
#                 small_batch = uids[uid_start:uid_end]
#                 if small_batch_forced is not None:
#                     small_batch = small_batch_forced
#                 if len(small_batch) == 0:
#                     break
#                 ##
#                 to_download_uids |= set(small_batch)
#                 ##
#                 small_batch_count = len(small_batch)
#                 small_batch_text = ','.join(small_batch)
#                 epost = Entrez.read(Entrez.epost(db, id=small_batch_text))
#                 webenv = epost['WebEnv']
#                 query_key = epost['QueryKey']

#                 temp_records = []

#                 for start in range(0, small_batch_count, small_batch_size):
#                     end = min(small_batch_count, start + small_batch_size)
#                     print ('  Going to download record %i to %i of %i.'
#                            % (start + 1, end, small_batch_count))

#                     # for i, j in enumerate(range(start, end)):
#                     #     print(i, small_batch[j])

#                     fetch_handle = Entrez.efetch(
#                         db=db, rettype=rettype, retmode=retmode,
#                         retstart=start, retmax=small_batch_size,
#                         webenv=webenv, query_key=query_key)

#                     batch_data = krbioio.read_sequence_data(fetch_handle, rettype)
#                     temp_records = temp_records + batch_data

#                 # n_rec_to_download = uid_end - uid_start
#                 n_rec_to_download = len(small_batch)
#                 rec_downloaded = len(temp_records)
#                 ##
#                 for tr in temp_records:
#                     downloaded_uids.add(tr.annotations['gi'])
#                 ##
#             except Exception as err:
#                 # print(rec_downloaded, n_rec_to_download)
#                 # print(len(downloaded_uids), len(to_download_uids))
#                 # print(downloaded_uids - to_download_uids)
#                 # print(to_download_uids - downloaded_uids)

#                 print(err)

#                 # print('    HTTP problem, retrying...')
#                 # time.sleep(5)
#                 # continue

#             # print(rec_downloaded, n_rec_to_download)
#             # print(len(downloaded_uids), len(to_download_uids))
#             # print(downloaded_uids - to_download_uids)
#             # print(to_download_uids - downloaded_uids)

#             if rec_downloaded == n_rec_to_download:
#                 print('    Downloaded', rec_downloaded, 'of',
#                       n_rec_to_download, 'records.')
#                 SeqIO.write(temp_records, out_handle, 'gb')
#                 fetch_handle.close()
#                 small_batch_forced = None
#                 break
#             else:
#                 fetch_handle.close()
#                 missing_uids_now = to_download_uids - downloaded_uids
#                 if len(missing_uids_now & missing_uids) > 0:
#                     print('    Excluding uids:', ', '.join(list(missing_uids_now)))
#                     small_batch_forced = list(set(small_batch) - missing_uids_now)
#                     continue
#                 # print('    Downloaded', rec_downloaded, 'of',
#                 #       n_rec_to_download, 'records.')
#                 print('    Download corrupted, retrying...')
#                 print('    Missing uids:', ', '.join(list(missing_uids_now)))
#                 missing_uids |= missing_uids_now
#                 continue

#     out_handle.close()

#     return


def get_ncbi_tax_id_for_record(record):
    import krseq
    feature_index = krseq.get_features_with_qualifier(
        record, 'db_xref', 'taxon', feature_type=None, loose=True)[0]
    # ToDo: search for taxon key as we now assume that it will always be the
    # last in the list
    return int(record.features[feature_index].qualifiers['db_xref'][-1].split('taxon:')[1])


def get_ncbi_tax_id_for_tax_term(email, tax_term):

    taxid_list = list(esearch(tax_term, 'taxonomy', email))
    taxid = None
    if len(taxid_list) > 0:
        taxid = int(taxid_list[0])

    return taxid


def get_taxids(email, tax_terms=None, tax_ids=None):

    # from Bio import Entrez

    if not tax_ids:
        tax_ids = list()
        for tax_term in tax_terms:
            tax_id = get_ncbi_tax_id_for_tax_term(email, tax_term)
            if tax_id:
                tax_ids.append(tax_id)

    # clean_tax_ids = list()
    # for ti in tax_ids:
    #     clean_tax_ids.append(str(ti))

    records = list()
    if tax_ids:
        # Entrez.email = email
        # handle = Entrez.efetch('taxonomy', id=','.join(clean_tax_ids), retmode="xml")
        # records = Entrez.read(handle)

        records = download_ncbi_records(
            uids=tax_ids,
            db='taxonomy',
            entrez_email=email,
            rettype=None,
            large_batch_size=2000,
            small_batch_size=500,
            genbank_mode=False)

    results = dict()
    for i, record in enumerate(records):
        ###
        # for kk in record.keys():
        #     print(kk, '  ::  ', record[kk])
        #     print('\n\n')
        ###
        # lineage_list_temp = record['LineageEx']
        # lineage_list = list()
        # for l in lineage_list_temp:
        #     lin_item = 'name=' + l['ScientificName'] + ';' + 'rank=' + l['Rank'] + ';' + 'taxid=' + l['TaxId'] + ';'
        #     lineage_list.append(lin_item)
        results_rec = list()
        tax_id_main = record['TaxId']
        results_rec.append(tax_id_main)
        if 'AkaTaxIds' in record.keys():
            other_tax_ids = record['AkaTaxIds']
            results_rec = results_rec + other_tax_ids
        results_rec = [int(x) for x in results_rec]
        if tax_terms:
            results[str(tax_terms[i])] = results_rec
        else:
            results[str(tax_ids[i])] = results_rec

    return results


def get_lineages(email, tax_terms=None, tax_ids=None):

    # from Bio import Entrez

    if not tax_ids:
        tax_ids = list()
        for tax_term in tax_terms:
            tax_id = get_ncbi_tax_id_for_tax_term(email, tax_term)
            if tax_id:
                tax_ids.append(tax_id)

    # clean_tax_ids = list()
    # for ti in tax_ids:
    #     clean_tax_ids.append(str(ti))

    records = list()
    if tax_ids:
        # Entrez.email = email
        # handle = Entrez.efetch('taxonomy', id=','.join(clean_tax_ids), retmode="xml")
        # records = Entrez.read(handle)

        records = download_ncbi_records(
            uids=tax_ids,
            db='taxonomy',
            entrez_email=email,
            rettype=None,
            large_batch_size=2000,
            small_batch_size=500,
            genbank_mode=False)

    results = dict()
    for record in records:
        ###
        # for kk in record.keys():
        #     print(kk, '  ::  ', record[kk])
        #     print('\n\n')
        ###
        lineage_list_temp = record['LineageEx']
        lineage_list = list()
        for l in lineage_list_temp:
            lin_item = 'name=' + l['ScientificName'] + ';' + 'rank=' + l['Rank'] + ';' + 'taxid=' + l['TaxId'] + ';'
            lineage_list.append(lin_item)
        tax_id = record['TaxId']
        results[tax_id] = lineage_list

    return results


def parse_lineage_string_list(lineage_string_list):

    full_list = list()
    name_list = list()

    for ls in lineage_string_list:
        ls_split_1 = ls.split(';')
        local_dict = dict()
        for ls1 in ls_split_1:
            ls_split_2 = ls1.split('=')
            if ls_split_2[0] == '':
                continue
            # print(ls_split_2)
            key = ls_split_2[0]
            value = ls_split_2[1]
            if key == 'name':
                name_list.append(value)
            local_dict[key] = value
        full_list.append(local_dict)
        #     print(key, value)
        # print('--- --- --- --- ---')

    return [full_list, name_list]


def get_common_names(email, tax_ids=None):

    # from Bio import Entrez

    # clean_tax_ids = list()
    # for ti in tax_ids:
    #     clean_tax_ids.append(str(ti))

    # Entrez.email = email
    # handle = Entrez.efetch('taxonomy', id=','.join(clean_tax_ids), retmode="xml")
    # records = Entrez.read(handle)

    records = download_ncbi_records(
        uids=tax_ids,
        db='taxonomy',
        entrez_email=email,
        rettype=None,
        large_batch_size=2000,
        small_batch_size=500,
        genbank_mode=False)

    results = dict()
    for record in records:
        for k in record.keys():
            # print(k, ' :: ', record[k])
            # print('--- --- ---')
            if ('OtherNames' in record.keys()) and ('GenbankCommonName' in record['OtherNames'].keys()):
                common_name = record['OtherNames']['GenbankCommonName']
                tax_id = record['TaxId']
                results[tax_id] = common_name
        # print('================================================================')

    return results


def get_common_name(email, tax_term=None, tax_id=None):

    # from Bio import Entrez

    taxid = None
    if not tax_id:
        taxid = get_ncbi_tax_id_for_tax_term(email, tax_term)
    else:
        taxid = tax_id

    common_name = None
    if taxid:
        # Entrez.email = email
        # handle = Entrez.efetch('taxonomy', id=str(taxid), retmode="xml")
        # record = Entrez.read(handle)[0]

        record = download_ncbi_records(
            uids=[taxid],
            db='taxonomy',
            entrez_email=email,
            rettype=None,
            large_batch_size=2000,
            small_batch_size=500,
            genbank_mode=False)[0]

        # for k in record['OtherNames'].keys():
        #     print(k, ' :: ', record['OtherNames'][k])
        #     print('--- --- ---')

        common_name = record['OtherNames']['GenbankCommonName']

    return common_name


# if __name__ == '__main__':

    # Tests

    # dbs = entrez_db_list(email='test@test.com')
    # print(dbs)

    # download_ncbi_records(
    #     uids=[379142464, 31281443, 31281445],
    #     db='nuccore',
    #     entrez_email='test@test.com',
    #     rettype='gb',
    #     large_batch_size=2000,
    #     small_batch_size=500,
    #     genbank_mode=False)

    # download_sequence_records(
    #     file_path='/Users/karolis/Desktop/test.gb',
    #     uids=[379142464, 31281443, 31281445],
    #     db='nuccore',
    #     entrez_email='test@test.com')

    # taxids = get_taxids(
    #     email='test@test.com',
    #     tax_terms=None,
    #     tax_ids=[9606])
    # print(taxids)

    # common_name = get_common_name(email='test@test.com', tax_term='Callithrix geoffroyi')
    # print(common_name)

    # common_names = get_common_names(email='test@test.com', tax_ids=['52231', '9483'])
    # print(common_names)

