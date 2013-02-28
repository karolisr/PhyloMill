from __future__ import print_function
#from __future__ import unicode_literals


def parse_organism_name(name, sep=' ', ncbi_authority=False):

    '''
    Parse an organism name and return a dictionary with fields:

        genus
        species
        authority
        variety
        subspecies
        cross
        other - eg.: Solanum sp. 112
    '''

    var_string = 'var.'
    form_string = 'f.'
    sub_string = 'subsp.'
    cross_string = 'x'
    other_string = 'sp.'

    name_list = name.split(sep)

    var_bool = bool(name_list.count(var_string))
    form_bool = bool(name_list.count(form_string))
    sub_bool = bool(name_list.count(sub_string))
    cross_bool = bool(name_list.count(cross_string))
    other_bool = bool(name_list.count(other_string))

    var = None
    form = None
    sub = None
    cross = None
    other = None

    if var_bool:
        var = name_list.index(var_string)
    if form_bool:
        form = name_list.index(form_string)
    if sub_bool:
        sub = name_list.index(sub_string)
    if cross_bool:
        cross = name_list.index(cross_string)
    if other_bool:
        other = name_list.index(other_string)

    var_dict = {'name': 'variety', 'index': var}
    form_dict = {'name': 'form', 'index': form}
    sub_dict = {'name': 'subspecies', 'index': sub}
    cross_dict = {'name': 'cross', 'index': cross}
    other_dict = {'name': 'other', 'index': other}

    indexes = []
    if var is not None:
        indexes.append(var_dict)
    if form is not None:
        indexes.append(form_dict)
    if sub is not None:
        indexes.append(sub_dict)
    if cross is not None:
        indexes.append(cross_dict)
    if other is not None:
        indexes.append(other_dict)

    indexes.sort(key=lambda x: x['index'])

    organism_dict = {}

    organism_dict['genus'] = name_list[0]
    organism_dict['species'] = ''
    organism_dict['authority'] = ''

    if ncbi_authority:
        if ((len(name_list) > 2) and
            (len(indexes) and
             indexes[0]['index'] >= 2)):
            organism_dict['authority'] = (
                ' '.join(name_list[indexes[-1]['index'] + 2:len(name_list)]))
        if len(name_list) > 2 and len(indexes) == 0:
            organism_dict['authority'] = ' '.join(name_list[2:len(name_list)])

    if (len(name_list) >= 2 and ((len(indexes) == 0) or
           (len(indexes) and indexes[0]['index'] != 1))):
        organism_dict['species'] = name_list[1]

    if len(indexes) and indexes[0]['index'] > 1:
        organism_dict['species'] = name_list[1]

    organism_dict['variety'] = ''
    organism_dict['form'] = ''
    organism_dict['subspecies'] = ''
    organism_dict['cross'] = ''
    organism_dict['other'] = ''

    number_of_indexes = len(indexes)

    if number_of_indexes:
        previous_index = None
        for i, index in enumerate(indexes):
            if previous_index is None:
                previous_index = index
                continue

            '''
            NCBI authority information is a mess. Sometimes it is appended to
            the binomial in the 'authority' 'name_class', however sometimes
            it is appended to the 'synonym' 'name_class'. If it follows a
            variety or subspecies designation there is no way to tell when this
            designation ends and the authority information begins. A heuristic
            I am using here is that var. and subsp. are usually a single word.
            There is a problem however when the species named is a cross 'x',
            e.g. 'Genus x specific' or 'Genus specific x Genus specific'.
            Sometimes x is followed by two words. Therefore, if we are
            expecting to see NCBI style appended authority information, we
            should set ncbi_authority=True. However, some cross species will
            be parsed incorrectly. If no NCBI authority information is
            expected, ncbi_authority=False (default) should be used. In case
            there is NCBI authority information after all, it will be appended
            to the preceding field.
            '''

            # Note: For our purposes crosses are bad anyways, right?!
            # So we don't much care.

            if ncbi_authority:
                organism_dict[previous_index['name']] = (
                    ' '.join(name_list[previous_index['index'] +
                                       1:previous_index['index'] + 2]))
            else:
                organism_dict[previous_index['name']] = (
                    ' '.join(name_list[previous_index['index'] +
                                       1:index['index']]))
            previous_index = index

        if ncbi_authority:
            organism_dict[previous_index['name']] = (
                ' '.join(name_list[previous_index['index'] +
                                   1:previous_index['index'] + 2]))
        else:
            organism_dict[previous_index['name']] = (
                ' '.join(name_list[previous_index['index'] +
                                   1:len(name_list)]))

    return organism_dict


def flatten_organism_name(parsed_name, sep=' '):

    '''
    Take a dictionary from "parse_organism_name" and return an organism name
    string.
    '''

    genus_bool = 'genus' in parsed_name
    species_bool = 'species' in parsed_name
    var_bool = 'variety' in parsed_name
    form_bool = 'form' in parsed_name
    sub_bool = 'subspecies' in parsed_name
    cross_bool = 'cross' in parsed_name
    other_bool = 'other' in parsed_name

    name = ''

    var_string = 'var.'
    form_string = 'f.'
    sub_string = 'subsp.'
    cross_string = 'x'
    other_string = 'sp.'

    if genus_bool and parsed_name['genus'] != '':
        name = name + parsed_name['genus']
    if species_bool and parsed_name['species'] != '':
        name = name + sep + parsed_name['species']
    if cross_bool and parsed_name['cross'] != '':
        name = name + sep + cross_string + sep + parsed_name['cross']
    if other_bool and parsed_name['other'] != '':
        name = name + sep + other_string + sep + parsed_name['other']
    if var_bool and parsed_name['variety'] != '':
        name = name + sep + var_string + sep + parsed_name['variety']
    if form_bool and parsed_name['form'] != '':
        name = name + sep + form_string + sep + parsed_name['form']
    if sub_bool and parsed_name['subspecies'] != '':
        name = name + sep + sub_string + sep + parsed_name['subspecies']

    return name


def accepted_name(name, synonymy_table, auth_file, sep=' ',
                  allow_loose_matching=True, ncbi_authority=False, level=1):
    import copy
    # Emma Goldberg's module to standardize authority
    from krtp.eg import stdauth
    # Make a copy of the synonymy_table object. We will be removing entries
    # from it, we do not want this to affect the original object.
    synonymy_table = copy.copy(synonymy_table)
    authority_alternates = stdauth.make_auth_dic(auth_file)
    # Organism name parsed
    o = None
    if isinstance(name, basestring):
        o = parse_organism_name(name, sep=sep, ncbi_authority=ncbi_authority)
    else:
        o = name
    # Take available authority information and translate it into an
    # accepted form.
    o['authority'] = stdauth.translate(o['authority'], authority_alternates)
    #print(level * '-', 'o', str(level), flatten_organism_name(o))
    # This will store the accepted name
    accepted = dict()
    accepted['genus'] = ''
    accepted['species'] = ''
    accepted['authority'] = ''
    accepted['subspecies'] = ''
    accepted['variety'] = ''
    accepted['status'] = ''
    accepted['id'] = ''

    matching_entries_strict = list()
    matching_entries_loose = list()
    matching_entries = list()

    # Find the entries in synonymy table that match our organism name.
    for s in synonymy_table:
        if (o['genus'] == s['Genus'] and
            o['species'] == s['Species'] and
            o['variety'] == s['Variety'] and
            o['subspecies'] == s['Subspecies'] and
                o['authority'] == s['Authority']):
            matching_entries_strict.append(s)

    if allow_loose_matching:
        for s in synonymy_table:
            if (o['genus'] == s['Genus'] and
                o['species'] == s['Species'] and
                (o['variety'] == '' or s['Variety'] == '' or
                    (o['variety'] == s['Variety'])) and
                (o['subspecies'] == '' or s['Subspecies'] == '' or
                    (o['subspecies'] == s['Subspecies'])) and
                (o['authority'] == '' or s['Authority'] == '' or
                    (o['authority'] == s['Authority']))):
                matching_entries_loose.append(s)

    matching_entries_strict = sorted(matching_entries_strict, key=lambda x: (
        x['Authority'],
        x['Variety'],
        x['Subspecies']))
    matching_entries_strict.reverse()
    matching_entries_loose = sorted(matching_entries_loose, key=lambda x: (
        x['Authority'],
        x['Variety'],
        x['Subspecies']))
    matching_entries_loose.reverse()
    matching_entries = matching_entries_strict + matching_entries_loose

    #print(' ' * level, '=== === === ===')
    #for e in matching_entries:
    #    print(' ' * level, e['Genus'], e['Species'], 'A', e['Authority'],
    #        'S', e['Subspecies'], 'V', e['Variety'])
    #print(' ' * level, '=== === === ===')

    for s in matching_entries:
        accepted['genus'] = s['AccGenus']
        accepted['species'] = s['AccSpecies']
        accepted['authority'] = s['AccAuthority']
        accepted['subspecies'] = s['AccSubspecies']
        accepted['variety'] = s['AccVariety']
        accepted['status'] = s['Status']
        accepted['id'] = s['AccID']
        # If the matching entry is a synonym, recurse into synonymy table
        # until an entry with a non-synonym status is reached.
        if (s['Status'].lower() == 'syn' or
                s['Status'].lower() == 'syn-alt'):
            if s in synonymy_table:
                synonymy_table.remove(s)
            return(accepted_name(accepted, synonymy_table, auth_file,
                                 sep=sep, level=level + 1))
        else:
            #print(level * '-', 'a', str(level),
            #    flatten_organism_name(accepted))
            return(accepted)


def names_for_ncbi_taxid(tax_id, ncbi_names_table, sorting='class'):

    '''
    Return all the names ("synonyms") associated with an NCBI taxid.
    '''

    names = list()
    #sci_name = None
    #authority_name = None
    for row in ncbi_names_table:
        if row['tax_id'] == str(tax_id):
            names.append(row)

    auth_names = list()
    syn_names = list()
    sci_names = list()

    for row in names:
        parsed = parse_organism_name(row['name_txt'],
                                     ncbi_authority=True)

        # NCBI names table includes common names and other weird things, we do
        # not want any of that.

        if row['name_class'] == 'scientific name':
            parsed['name_class'] = '2'
            sci_names.append(parsed)
        if row['name_class'] == 'authority':
            parsed['name_class'] = '1'
            auth_names.append(parsed)
        if row['name_class'] == 'synonym':
            parsed['name_class'] = '3'
            syn_names.append(parsed)

    priority_list = auth_names + syn_names + sci_names
    # This will sort the names so the results with authority information (if
    # any) will appear at the beginning of the list.
    if sorting.startswith('class'):
        priority_list.sort(key=lambda x: x['name_class'], reverse=False)
    elif sorting.startswith('authority'):
        priority_list.sort(key=lambda x: x['authority'], reverse=True)
    return priority_list

if __name__ == '__main__':

    # Tests

    import os

    ps = os.path.sep

    # parse_organism_name
    #name = parse_organism_name('A b x c var. d subsp. e', sep=' ',
    #    ncbi_authority=False)
    #print(name)

    # flatten_organism_name
    #print(flatten_organism_name(name))

    # accepted_name
    #import krio
    #synonymy_table = krio.read_table_file('testdata' + ps + 'synonymy.csv',
    #    has_headers=True, headers=None, delimiter=',')

    #synonymy_table = krio.read_table_file(
    #    '/home/karolis/Dropbox/code/test/sol-in-1/synonymy',
    #    has_headers=True, headers=None, delimiter=',')

    #an1 = accepted_name('Solanum x juzepczukii', synonymy_table,
    #    'testdata' + ps + 'authority_alternates.dat', sep=' ')

    #an2 = accepted_name('Petunia axillaris (Lam.) Britton, Stern & Poggenb.',
    #    synonymy_table, 'testdata' + ps + 'authority_alternates.dat', sep=' ')

    #Pprint(an1)
    #print(an2)

    # names_for_ncbi_taxid
    #ncbi_names_table = krio.read_table_file(
    #    'testdata' + ps + 'ncbi_names.dmp',
    #    has_headers=False,
    #    headers=('tax_id', 'name_txt', 'unique_name', 'name_class'),
    #    delimiter='\t|')
    #ncbi_names_table = krio.read_table_file(
    #    '/home/karolis/Dropbox/code/test/sol-in-1/ncbi-names',
    #    has_headers=False,
    #    headers=('tax_id', 'name_txt', 'unique_name', 'name_class'),
    #    delimiter='\t|')
    #ncbi_names_table = krio.read_table_file(
    #    path='/home/karolis/Dropbox/code/test/sol-in-1/ncbi-names',
    #    has_headers=False,
    #    headers=('tax_id', 'name_txt', 'unique_name',
    #        'name_class'),
    #    delimiter='\t|',
    #    quotechar=None,
    #    stripchar='"',
    #    rettype='dict')
    #names = names_for_ncbi_taxid('136614', ncbi_names_table, sorting='class')
    #for n in names:
    #    print(n)
