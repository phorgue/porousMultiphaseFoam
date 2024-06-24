#!/usr/bin/python

import yaml

yaml_file = open('../etc/keywords.yaml', 'r')
dict_file = open('../libraries/keywords/keywords_map.H', 'w')

data = yaml.load(yaml_file, Loader=yaml.FullLoader)

for i in data:
    for j in data[i]:
        if isinstance(data[i][j], dict):
            for k in data[i][j]:
                if isinstance(data[i][j][k], dict):
                    for l in data[i][j][k]:
                        dict_file.write('{"'+i+'_'+j+'_'+k+'_'+l+'", "'+data[i][j][k][l]+'"},\n')
                else:
                    dict_file.write('{"'+i+'_'+j+'_'+k+'", "'+data[i][j][k]+'"},\n')
                    
        else:
            dict_file.write('{"'+i+'_'+j+'", "'+data[i][j]+'"},\n')