chr22_dinuc_freq = {'AA': 0.07946777262528104,
                    'AC': 0.05093424568835639,
                    'AG': 0.07486997248981177,
                    'AT': 0.06053258477333608,
                    'CA': 0.07631952805132396,
                    'CC': 0.06695292775053881,
                    'CG': 0.016281693225376138,
                    'CT': 0.07386697509406903,
                    'GA': 0.06275311378792536,
                    'GC': 0.05427934466020802,
                    'GG': 0.06817321590127282,
                    'GT': 0.05092521567010435,
                    'TA': 0.047264446270726033,
                    'TC': 0.061224664382737416,
                    'TG': 0.07680695893128689,
                    'TT': 0.0793473406976459}

read_dinuc_freq = {'AA': 0.09177332285779456,
                   'AC': 0.04823884228400537,
                   'AG': 0.05506703500012672,
                   'AT': 0.08062599792178828,
                   'CA': 0.06801682844615657,
                   'CC': 0.07035963200446055,
                   'CG': 0.06304027168816687,
                   'CT': 0.052398307017766176,
                   'GA': 0.08243353524089515,
                   'GC': 0.061652431761157714,
                   'GG': 0.08126669538991814,
                   'GT': 0.03717363205514864,
                   'TA': 0.03282510074258053,
                   'TC': 0.07245811896494919,
                   'TG': 0.0643449831462098,
                   'TT': 0.038325265478875735}

read_nuc_freq = {'A': 0.27542949277043843,
                 'C': 0.2533452225680747,
                 'G': 0.2629793035295624,
                 'T': 0.2082459811319245}

chr22_nuc_freq = {'A': 0.26580368810109145,
                  'C': 0.23340643550863757,
                  'G': 0.23611855838025542,
                  'T': 0.2646713180100156}



rel_dinuc_freq_dict = {}

#for key in chr22_dinuc_freq:
#    rel_dinuc_freq_dict.setdefault('read_{} / chr22_{}'.format(key,key), read_dinuc_freq[key] / chr22_dinuc_freq[key])

#for key, value in rel_freq_dict.items():
#    print('{} : {}'.format(key, value))

rel_nuc_dict = {}

for key in read_nuc_freq:
    rel_nuc_dict.setdefault('read_{} / chr22_{}'.format(key,key), read_nuc_freq[key] / chr22_nuc_freq[key])

for key, value in rel_nuc_dict.items():
    print('{} : {}'.format(key, value))
