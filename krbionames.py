from __future__ import print_function
from __future__ import unicode_literals

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
    if var != None:
        indexes.append(var_dict)
    if form != None:
        indexes.append(form_dict)
    if sub != None:
        indexes.append(sub_dict)
    if cross != None:
        indexes.append(cross_dict)
    if other != None:
        indexes.append(other_dict)

    indexes.sort(key=lambda x:x['index'])

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
            if previous_index == None:
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
                                       1:previous_index['index']+2]))
            else:
                organism_dict[previous_index['name']] = (
                    ' '.join(name_list[previous_index['index'] +
                                       1:index['index']]))
            previous_index = index

        if ncbi_authority:
            organism_dict[previous_index['name']] = (
                ' '.join(name_list[previous_index['index'] +
                                   1:previous_index['index']+2]))
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
    
    genus_bool = parsed_name.has_key('genus')
    species_bool = parsed_name.has_key('species')
    var_bool = parsed_name.has_key('variety')
    form_bool = parsed_name.has_key('form')
    sub_bool = parsed_name.has_key('subspecies')
    cross_bool = parsed_name.has_key('cross')
    other_bool = parsed_name.has_key('other')

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

def accepted_name(name, synonymy_table, auth_file, sep=' '):

    '''
    Takes the organism name, either as a string or an output of
    "parse_organism_name" function and returns an accepted name based on
    synonymy_table information.
    '''

    # Emma Goldberg's module to standardize authority
    from krtp.eg import stdauth
    authority_alternates = stdauth.make_auth_dic(auth_file)

    # Organism name parsed
    o = None
    if isinstance(name, basestring):
        o = parse_organism_name(name)
    else:
        o = name
    # Take available authority information and translate it into an
    # accepted form.
    o['authority'] = stdauth.translate(o['authority'], authority_alternates)
    accepted = dict()
    accepted['genus'] = ''
    accepted['species'] = ''
    accepted['variety'] = ''
    accepted['subspecies'] = ''
    accepted['status'] = ''
    accepted['authority'] = ''
    accepted['id']= ''
    matching_entries = list()

    # Find the entries in synonymy table that match our organism name.
    for s in synonymy_table:
        if (o['genus'] == s['Genus'] and
            o['species'] == s['Species'] and
            o['variety'] == s['Variety'] and
            o['subspecies'] == s['Subspecies'] and
            (o['authority'] == '' or (o['authority'] == s['Authority']))
            ):
            matching_entries.append(s)

    for s in matching_entries:
        accepted['genus'] = s['AccGenus']
        accepted['species'] = s['AccSpecies']
        accepted['variety'] = s['AccVariety']
        accepted['subspecies'] = s['AccSubspecies']
        accepted['status'] = s['Status']
        accepted['authority'] = s['AccAuthority']
        accepted['id'] = s['AccID']
        # If the matching entry is a synonym, recurse into synonymy table
        # until an entry with a non-synonym status is reached.
        if (s['Status'].lower() == 'syn' or
            s['Status'].lower() == 'syn-alt'):
            accepted = accepted_name(accepted, synonymy_table, auth_file,
                                     sep=sep)
        else:
            break
    return accepted

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep

    # parse_organism_name
    name = parse_organism_name('A b x c var. d subsp. e', sep=' ',
        ncbi_authority=False)
    print(name)

    # flatten_organism_name
    print(flatten_organism_name(name))

    # accepted_name
    import krio
    synonymy_table = krio.read_table_file('testdata'+PS+'synonymy.csv',
        has_headers=True, headers=None, delimiter=b',', iterator=False)
    an = accepted_name('Physalis microphysa', synonymy_table,
        'testdata'+PS+'authorityalternates.dat', sep=b' ')
    print(an)