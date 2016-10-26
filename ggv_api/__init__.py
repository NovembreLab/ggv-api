from flask import Flask,json,Response
from flask.ext.restful import Api,Resource,reqparse
from flask.ext.restful.utils import cors
import subprocess
import os
import random
import httplib2
import sys
import json
from pyliftover import LiftOver

'''
A simple API to extract Fst and Freq info from tabix-accessible files
John Novembre and Joe Marcus
'''

# Basic header info
app=Flask(__name__, static_url_path='')
api=Api(app)
api.decorators=[cors.crossdomain(origin='*')]

# Basic app route
@app.route('/')
def index():
    return "API for accessing population frequency data from Novembre Lab."

# Define parser for arguments
parser = reqparse.RequestParser()
parser.add_argument('chr',type=int)
parser.add_argument('pos',type=int)
parser.add_argument('rsID',type=str)
parser.add_argument('random_snp',type=bool)
parser.add_argument('data',type=str,required=True)

# Converting rsID
lo = LiftOver('hg19', 'hg18')

# The data
data_files = {
                '"1000genomes_phase3_table"':{
                                  'data_dir':'/var/www/ggv_api/ggv_api/static/data/1000genomes_phase3_table/',
                                  'file_base':'ALL.chr',
                                  'file_end':'.phase3_shapeit2_mvncall_integrated.20130502.genotype.freq.gz',
                                  'random_end':'.phase3_shapeit2_mvncall_integrated.20130502.genotype.snps.txt',
                                  'pop_file':'1000genomes_phase3_clst.txt',
                                  'coordinates':'1000genomes_phase3_coordinates.txt',
                                  'multiple_chr':True
                                 },
                '"hgdp_table"':{
                                 'data_dir':'/var/www/ggv_api/ggv_api/static/data/H938_table/',
                                 'file':'H938_autoSNPs_full.freq.gz',
                                 'random_file':'H938_autoSNPs_full.snps.txt',
                                 'pop_file':'H938_pops_full.txt',
                                 'coordinates':'H938_coordinates.txt',
                                 'multiple_chr':False
                                },
                '"ExAC_table"':{
                                 'data_dir':'/var/www/ggv_api/ggv_api/static/data/ExAc_table/',
                                 'file_base':'ExAC_chr',
                                 'file_end':'.frq.gz',
                                 'random_end':'.snps.txt',
                                 'pop_file':'ExAc_pops.txt',
                                 'coordinates':'ExAc_coordinates.txt',
                                 'multiple_chr':True
                                },
                '"1KG_phase3_superpops_table"':{
                                 'data_dir':'/var/www/ggv_api/ggv_api/static/data/1KG_phase3_superpops_table/',
                                 'file_base':'1KG_superpop_chr',
                                 'file_end':'.frq.gz',
                                 'random_end':'.snps.txt',
                                 'pop_file':'1KG_superpop_pops.txt',
                                 'coordinates':'1KG_superpop_coordinates.txt',
                                 'multiple_chr':True
                                },
                '"popres_euro_table"':{
                                  'data_dir':'/var/www/ggv_api/ggv_api/static/data/POPRES_EURO_hg19_table/',
                                  'file_base':'POPRES_NovembreEtAl2008_autoSNPs_hg19.freq_chr',
                                  'file_end':'.srt.gz',
                                  'random_end':'.snps.txt',
                                  'pop_file':'POPRES_euro_pops.txt',
                                  'coordinates':'POPRES_euro_coordinates.txt',
                                  'multiple_chr':True
                                 }
             }


class FreqTable(object):
    '''
    class for obtaining frequencies from FreqVcf Tables
    '''
    def __init__(self,dataset):
        self.dataset=dataset

    def get_by_chr_pos(self, chr, pos):
        '''
        '''
        if data_files[self.dataset]['multiple_chr']:
            vcf_filename = data_files[self.dataset]['data_dir']+data_files[self.dataset]['file_base']+str(chr)+data_files[self.dataset]['file_end']
        else:
            vcf_filename = data_files[self.dataset]['data_dir']+data_files[self.dataset]['file']

        tabix_snp_pos = str(chr)+':'+str(pos)+'-'+str(pos)
        tabix_command = '/home/josephhmarcus/bin/tabix-0.2.6/tabix '+vcf_filename+' '+tabix_snp_pos
        proc = subprocess.Popen(tabix_command, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        freq_dict = {}
        freq_list = out.strip('\n').split('\n')
        for line in freq_list:
            line = line.split('\t')
            freq_dict[line[3]] = line

        return self.freq_out_to_json(freq_dict)

    def get_by_rsID(self, rsID):
        '''
        '''
        rest_command = "curl 'http://grch37.rest.ensembl.org/variation/human/"+rsID+"?' -H 'Content-type:application/json'"
        proc = subprocess.Popen(rest_command, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        rsID_data = json.loads(out)
        snp_pos = str(rsID_data['mappings'][0]['location'])
        chr = snp_pos.split(':')[0]
        if data_files[self.dataset]['multiple_chr']:
            tabix_snp_pos = snp_pos
            vcf_filename = data_files[self.dataset]['data_dir']+data_files[self.dataset]['file_base']+str(chr)+data_files[self.dataset]['file_end']
        else:
            pos = snp_pos.split('-')[1]
            chrom = 'chr'+chr
            lo_list = lo.convert_coordinate(chrom, int(pos))
            tabix_snp_pos = lo_list[0][0][3:]+':'+str(lo_list[0][1])+'-'+str(lo_list[0][1])
            vcf_filename = data_files[self.dataset]['data_dir']+data_files[self.dataset]['file']

        tabix_command = '/home/josephhmarcus/bin/tabix-0.2.6/tabix '+vcf_filename+' '+tabix_snp_pos
        proc = subprocess.Popen(tabix_command, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        freq_dict = {}
        freq_list = out.strip('\n').split('\n')
        for line in freq_list:
            line = line.split('\t')
            freq_dict[line[3]] = line

        return self.freq_out_to_json(freq_dict)

    def get_random(self):
        '''
        '''
        if data_files[self.dataset]['multiple_chr']:
            chr = random.randint(1,22)
            snp_filename = data_files[self.dataset]['data_dir']+data_files[self.dataset]['file_base']+str(chr)+data_files[self.dataset]['random_end']
        else:
            snp_filename = data_files[self.dataset]['data_dir']+data_files[self.dataset]['random_file']
        snp_proc = subprocess.Popen('shuf -n 1 '+snp_filename, stdout=subprocess.PIPE, shell=True)
        (snp_out, snp_err) = snp_proc.communicate()
        arg_list = snp_out.strip('\n').split(' ')
        tabix_snp_pos = arg_list[0]+':'+arg_list[1]+'-'+arg_list[1]

        if data_files[self.dataset]['multiple_chr']:
            vcf_filename = data_files[self.dataset]['data_dir']+data_files[self.dataset]['file_base']+str(chr)+data_files[self.dataset]['file_end']
        else:
            vcf_filename = data_files[self.dataset]['data_dir']+data_files[self.dataset]['file']

        tabix_command = '/home/josephhmarcus/bin/tabix-0.2.6/tabix '+vcf_filename+' '+tabix_snp_pos
        proc = subprocess.Popen(tabix_command, stdout=subprocess.PIPE, shell=True)
        freq_dict = {}
        (out, err) = proc.communicate()
        freq_list = out.strip('\n').split('\n')
        for line in freq_list:
            line = line.split('\t')
            #sys.stderr.write(str(line))
            freq_dict[line[3]] = line

        return self.freq_out_to_json(freq_dict)

    def get_lon_lat_dict(self):
        '''
        helper function to read coords file and return dict of
        lat lons...
        '''
        with open(data_files[self.dataset]['data_dir']+data_files[self.dataset]['coordinates'], 'r') as coordinate_file:
            lon_lat_dict = {}
            coordinate_file.next()
            for line in coordinate_file:
                line = line.strip('\n').split(' ')
                lon_lat_dict[line[0]] = [line[1], line[2]]

        return lon_lat_dict

    def define_freqscale(self, json_data):
        '''
        helper function to add freqscale to json
        '''
        maxfreq=0
        for i in range(0,len(json_data)):
            if json_data[i]['rawfreq']>maxfreq:
                maxfreq=json_data[i]['rawfreq']
        if maxfreq<0.001:
            for i in range(0,len(json_data)):
                json_data[i]['freq']=[json_data[i]['rawfreq']*1000,1-json_data[i]['rawfreq']*1000]
                json_data[i]['freqscale']=0.001
        elif maxfreq<0.01:
            for i in range(0,len(json_data)):
                json_data[i]['freq']=[json_data[i]['rawfreq']*100,1-json_data[i]['rawfreq']*100]
                json_data[i]['freqscale']=0.01
        elif maxfreq<0.1:
            for i in range(0,len(json_data)):
                json_data[i]['freq']=[json_data[i]['rawfreq']*10,1-json_data[i]['rawfreq']*10]
                json_data[i]['freqscale']=0.1
        else:
            for i in range(0,len(json_data)):
                json_data[i]['freq']=[json_data[i]['rawfreq'],1-json_data[i]['rawfreq']]
                json_data[i]['freqscale']=1

        return json_data

    def freq_out_to_json(self, freq_dict):
        '''
        writes json data
        '''
        json_data =[]
        lon_lat_dict = self.get_lon_lat_dict()
        for pop in freq_dict.keys():
            map_pos = lon_lat_dict[pop]
            nobs =  freq_dict[pop][6]
            xobs = freq_dict[pop][7]
            freq = freq_dict[pop][8]
            chr_pos = str(freq_dict[pop][0])+':'+str(freq_dict[pop][1])
            alleles = [freq_dict[pop][4], freq_dict[pop][5]]
            if int(nobs) == 0:
                nobs = 'M'
                xobs = 'M'
            json_data.append({'pop':pop, 'pos':map_pos, 'nobs':nobs,
                              'rawfreq':float(freq), 'chrom_pos':chr_pos, 'alleles':alleles,
                              'xobs':xobs})

            json_data = self.define_freqscale(json_data)

        return json_data

class FreqVcfTable(Resource):
    def get(self):
        args=parser.parse_args()
        f = FreqTable(args['data'])
        if args['chr'] and args['pos']:
            return f.get_by_chr_pos(args['chr'], args['pos'])
        elif args['random_snp'] == True:
            return f.get_random()
        elif args['rsID']:
            return f.get_by_rsID(args['rsID'])
        else:
            return 403

# Set up API resources
api.add_resource(FreqVcfTable,'/freq_table')

if __name__ == '__main__':
    app.run(debug = True)
