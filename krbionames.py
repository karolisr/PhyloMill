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
    sub_string = 'subsp.'
    cross_string = 'x'
    other_string = 'sp.'

    name_list = name.split(sep)

    var_bool = bool(name_list.count(var_string))
    sub_bool = bool(name_list.count(sub_string))
    cross_bool = bool(name_list.count(cross_string))
    other_bool = bool(name_list.count(other_string))

    var = None
    sub = None
    cross = None
    other = None

    if var_bool:
        var = name_list.index(var_string)
    if sub_bool:
        sub = name_list.index(sub_string)
    if cross_bool:
        cross = name_list.index(cross_string)
    if other_bool:
        other = name_list.index(other_string)

    var_dict = {'name': 'variety', 'index': var}
    sub_dict = {'name': 'subspecies', 'index': sub}
    cross_dict = {'name': 'cross', 'index': cross}
    other_dict = {'name': 'other', 'index': other}

    indexes = []
    if var != None:
        indexes.append(var_dict)
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

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep

    # parse_organism_name
    print(parse_organism_name('A b x c var. d subsp. e',
        sep=' ', ncbi_authority=False))
