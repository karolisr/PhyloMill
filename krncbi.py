from __future__ import print_function
from __future__ import unicode_literals

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
            http://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options

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

    from Bio import Entrez
    Entrez.email = email
    if isinstance(esearch_terms, basestring):
        esearch_terms = [esearch_terms]
    retmax = None
    uid_set = []
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
    return uid_set

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep

    # entrez_db_list
    print(entrez_db_list('test@test.com'))

    # esearch
    print(esearch('GBSSI[Gene Name] AND txid4070[Organism]', 'nuccore',
        'test@test.com'))
