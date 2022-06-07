from run_tree_builder_api2 import *
params = json.load(open('params.json','r'))
print ('PARAMS:\n', params)
boutiques = pd.read_csv('../data/boutique+named.csv').loc[:1000, 'smiles'].values
prioritizers=['bkms,reaxys','bkms', 'reaxys']
save_prefix='6Jun2022_boutique_1000'

host_ip = ['35.238.247.50', '35.222.173.30', '34.136.208.25', '35.202.110.39', '34.133.210.74']
hosts = ['https://{}/'.format(ip) for ip in host_ip]

submit_parallel_job(boutiques, prioritizers, hosts, params, save_prefix, return_results=False)

