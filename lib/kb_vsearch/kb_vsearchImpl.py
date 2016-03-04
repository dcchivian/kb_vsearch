#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat
import numpy as np
import gzip

from Bio import SeqIO
from biokbase.workspace.client import Workspace as workspaceService
from requests_toolbelt import MultipartEncoder  # added
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService  # added

#END_HEADER


class kb_vsearch:
    '''
    Module Name:
    kb_vsearch

    Module Description:
    ** A KBase module: kb_vsearch
**
** This module contains 4 methods from VSEARCH: basic query/db search, clustering, chimera detection, and dereplication.
** 
** Initially only basic query/db search will be implemented between read sets
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None

    VSEARCH = '/kb/module/vsearch/bin/vsearch'

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def get_single_end_read_library(self, ws_data, ws_info, forward):
        pass

    def get_feature_set_seqs(self, ws_data, ws_info):
        pass

    def get_genome_feature_seqs(self, ws_data, ws_info):
        pass

    def get_genome_set_feature_seqs(self, ws_data, ws_info):
        pass

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


    # Helper script borrowed from the transform service, logger removed
    #
    def upload_file_to_shock(self,
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """

        if token is None:
            raise Exception("Authentication token required!")

        #build the header
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)
        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        #logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
        try:
            response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise
        if not response.ok:
            response.raise_for_status()
        result = response.json()
        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]


    def upload_SingleEndLibrary_to_shock_and_ws (self,
                                                 ctx,
                                                 workspace_name,
                                                 obj_name,
                                                 file_path,
                                                 provenance):

        # 1) upload files to shock
        token = ctx['token']
        forward_shock_file = self.upload_file_to_shock(
            shock_service_url = self.shockURL,
            filePath = file_path,
            token = token
            )
        #pprint(forward_shock_file)

        # 2) create handle
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        
        # 3) save to WS
        single_end_library = {
            'lib': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fasta',
                'size':forward_shock_file['file']['size']
            },
        }
        ws = workspaceService(self.workspaceURL, token=ctx['token'])
        new_obj_info = ws.save_objects({
                        'workspace':workspace_name,
                        'objects':[
                            {
                                'type':'KBaseAssembly.SingleEndLibrary',
                                'data':single_end_library,
                                'name':obj_name,
                                'meta':{},
                                'provenance':provenance
                            }]
                        })
        return new_obj_info[0]

    #END_CLASS_HEADER



    def VSearch_BasicSearch(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN VSearch_BasicSearch
        console = []
        self.log(console,'Running VSearch_BasicSearch with params=')
        self.log(console, "\n"+pformat(params))
        report = 'Running VSearch_BasicSearch with params='
        report += "\n"+pformat(params)


        #### do some basic checks
        objref = ''
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'input_one_name' not in params:
            raise ValueError('input_one_name parameter is required')
        if 'input_many_name' not in params:
            raise ValueError('input_many_name parameter is required')
        if 'output_filtered_name' not in params:
            raise ValueError('output_filtered_name parameter is required')

        #### Get the input_one object
        ##
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_one_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            # Object Info Contents
            # absolute ref = info[6] + '/' + info[0] + '/' + info[4]
            # 0 - obj_id objid
            # 1 - obj_name name
            # 2 - type_string type
            # 3 - timestamp save_date
            # 4 - int version
            # 5 - username saved_by
            # 6 - ws_id wsid
            # 7 - ws_name workspace
            # 8 - string chsum
            # 9 - int size 
            # 10 - usermeta meta
            one_type_name = info[2].split('.')[1].split('-')[0]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_one can be SingleEndLibrary or FeatureSet)
        #
        #  Note: currently only support SingleEndLibrary
        #
        if one_type_name == 'SingleEndLibrary':
            try:
                if 'lib' in data:
                    one_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    one_forward_reads = data['handle']
                else:
                    self.log(console,"bad structure for 'one_forward_reads'")
                    raise ValueError("bad structure for 'one_forward_reads'")
                #if 'lib2' in data:
                #    reverse_reads = data['lib2']['file']
                #elif 'handle_2' in data:
                #    reverse_reads = data['handle_2']
                #else:
                #    reverse_reads={}

                ### NOTE: this section is what could be replaced by the transform services
                one_forward_reads_file_path = os.path.join(self.scratch,one_forward_reads['file_name'])
                one_forward_reads_file_handle = open(one_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(one_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(one_forward_reads['url']+'/node/'+one_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    one_forward_reads_file_handle.write(chunk)
                one_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE

                #if 'interleaved' in data and data['interleaved']:
                #    self.log(console, 'extracting forward/reverse reads into separate files')
                #    if re.search('gz', forward_reads['file_name'], re.I):
                #        bcmdstring = 'gunzip -c ' + forward_reads_file_path
                #    else:    
                #        bcmdstring = 'cat ' + forward_reads_file_path 
                # 
                #    cmdstring = bcmdstring + '| (paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" > '+self.scratch+'/forward.fastq) | cut -f 5-8 | tr "\t" "\n" > '+self.scratch+'/reverse.fastq )'
                #    cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
                #    stdout, stderr = cmdProcess.communicate()
                #
                #    self.log(console, "cmdstring: " + cmdstring + " stdout: " + stdout + " stderr: " + stderr)
                #    
                #    forward_reads['file_name']='forward.fastq'
                #    reverse_reads['file_name']='reverse.fastq'
                #else:
                #    ### NOTE: this section is what could also be replaced by the transform services
                #    reverse_reads_file_path = os.path.join(self.scratch,reverse_reads['file_name'])
                #    reverse_reads_file_handle = open(reverse_reads_file_path, 'w', 0)
                #    self.log(console, 'downloading reverse reads file: '+str(reverse_reads_file_path))
                #    r = requests.get(reverse_reads['url']+'/node/'+reverse_reads['id']+'?download', stream=True, headers=headers)
                #    for chunk in r.iter_content(1024):
                #        reverse_reads_file_handle.write(chunk)
                #    reverse_reads_file_handle.close()
                #    self.log(console, 'done')
                #    ### END NOTE
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        #elif one_type_name == 'FeatureSet':
        #    # retrieve sequences for features

        else:
            raise ValueError('Cannot yet handle input_one type of: '+type_name)            


        #### Get the input_many object
        ##
        many_forward_reads_file_compression = None
        #sequencing_tech = 'artificial reads'
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['input_many_name']}])
            data = objects[0]['data']
            info = objects[0]['info']
            many_type_name = info[2].split('.')[1].split('-')[0]
            if data['lib']['file']['file_name'][-3:] == ".gz":
                many_forward_reads_file_compression = 'gz'
            #if 'sequencing_tech' in data:
            #    sequencing_tech = data['sequencing_tech']

        except Exception as e:
            raise ValueError('Unable to fetch input_many_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # Handle overloading (input_many can be SingleEndLibrary, FeatureSet, Genome, or GenomeSet)
        #
        #  Note: currently only support SingleEndLibrary
        #
        if many_type_name == 'SingleEndLibrary':

            # DEBUG
            #for k in data:
            #    self.log(console,"SingleEndLibrary ["+k+"]: "+str(data[k]))

            try:
                if 'lib' in data:
                    many_forward_reads = data['lib']['file']
                elif 'handle' in data:
                    many_forward_reads = data['handle']
                else:
                    self.log(console,"bad structure for 'many_forward_reads'")
                    raise ValueError("bad structure for 'many_forward_reads'")
                #if 'lib2' in data:
                #    reverse_reads = data['lib2']['file']
                #elif 'handle_2' in data:
                #    reverse_reads = data['handle_2']
                #else:
                #    reverse_reads={}

                ### NOTE: this section is what could be replaced by the transform services
                many_forward_reads_file_path = os.path.join(self.scratch,many_forward_reads['file_name'])
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'w', 0)
                self.log(console, 'downloading reads file: '+str(many_forward_reads_file_path))
                headers = {'Authorization': 'OAuth '+ctx['token']}
                r = requests.get(many_forward_reads['url']+'/node/'+many_forward_reads['id']+'?download', stream=True, headers=headers)
                for chunk in r.iter_content(1024):
                    many_forward_reads_file_handle.write(chunk)
                many_forward_reads_file_handle.close();
                self.log(console, 'done')
                ### END NOTE
            except Exception as e:
                print(traceback.format_exc())
                raise ValueError('Unable to download single-end read library files: ' + str(e))

        #elif many_type_name == 'Genome':
        #    # retrieve sequences for features

        #elif many_type_name == 'GenomeSet':
        #    # retrieve sequences for features

        else:
            raise ValueError('Cannot yet handle input_many type of: '+type_name)            


        ### Construct the command
        #
        #  e.g. vsearch --usearch_global data/input_many.fna --db data/input_one.fna --alnout output/alnout.txt --maxaccepts 10000 --maxrejects 100000000000 --wordlength 8 --minwordmatches 10 --id 0.5 -iddef 2
        #
        vsearch_cmd = [self.VSEARCH]

        # check for necessary files
        if not os.path.isfile(self.VSEARCH):
            raise ValueError("no such file '"+self.VSEARCH+"'")
        if not os.path.isfile(one_forward_reads_file_path):
            raise ValueError("no such file '"+one_forward_reads_file_path+"'")
        if not os.path.isfile(many_forward_reads_file_path):
            raise ValueError("no such file '"+many_forward_reads_file_path+"'")

        # set the output path
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch,'output.'+str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_aln_file_path = os.path.join(output_dir, 'alnout.txt');
        output_filtered_fasta_file_path = os.path.join(output_dir, 'output_filtered.fna');

        # this is command for basic search mode
        vsearch_cmd.append('--usearch_global')
        vsearch_cmd.append(many_forward_reads_file_path)
        vsearch_cmd.append('--db')
        vsearch_cmd.append(one_forward_reads_file_path)
        vsearch_cmd.append('--alnout')
        vsearch_cmd.append(output_aln_file_path)

        # options
        if 'maxaccepts' in params:
            if params['maxaccepts']:
                vsearch_cmd.append('--maxaccepts')
                vsearch_cmd.append(str(params['maxaccepts']))
        if 'maxrejects' in params:
            if params['maxrejects']:
                vsearch_cmd.append('--maxrejects')
                vsearch_cmd.append(str(params['maxrejects']))
        if 'wordlength' in params:
            if params['wordlength']:
                vsearch_cmd.append('--wordlength')
                vsearch_cmd.append(str(params['wordlength']))
        if 'minwordmatches' in params:
            if params['minwordmatches']:
                vsearch_cmd.append('--minwordmatches')
                vsearch_cmd.append(str(params['minwordmatches']))
        if 'ident_thresh' in params:
            if params['ident_thresh']:
                vsearch_cmd.append('--id')
                vsearch_cmd.append(str(params['ident_thresh']))
        if 'ident_mode' in params:
            if params['ident_mode']:
                vsearch_cmd.append('--iddef')
                vsearch_cmd.append(str(params['ident_mode']))


        # Run VSEARCH, capture output as it happens
        #
        self.log(console, 'running vsearch:')
        self.log(console, '    '+' '.join(vsearch_cmd))
        report += "\n"+'running vsearch:'+"\n"
        report += '    '+' '.join(vsearch_cmd)+"\n"

        p = subprocess.Popen(vsearch_cmd, \
                             cwd = self.scratch, \
                             stdout = subprocess.PIPE, \
                             stderr = subprocess.STDOUT, \
                             shell = False)

        while True:
            line = p.stdout.readline()
            if not line: break
            self.log(console, line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running vsearch, return code: '+str(p.returncode) + 
                '\n\n'+ '\n'.join(console))


        # Parse the valign output and store ids to filter many set to make filtered object to save back to KBase
        #
        # Valign calls hits "Query" (as strange as that is)
        #
        self.log(console, 'parsing vsearch alignment output')
        hit_seq_ids = dict()
        output_aln_file_handle = open (output_aln_file_path, "r", 0)
        output_aln_buf = output_aln_file_handle.readlines()
        output_aln_file_handle.close()
        hit_total = 0
        for line in output_aln_buf:
            # hits have lines of format 'Query >1367929'
            if line.startswith('Query >'):
                #self.log(console,'HIT LINE: '+line)  # DEBUG
                hit_total += 1
                hit_seq_id = line[7:]  # removes leading '>'
                if "\n" in hit_seq_id:
                    hit_seq_id = hit_seq_id[0:hit_seq_id.find("\n")+1]
                if "\t" in hit_seq_id:
                    hit_seq_id = hit_seq_id[0:hit_seq_id.find("\t")+1]
                if " " in hit_seq_id:
                    hit_seq_id = hit_seq_id[0:hit_seq_id.find(" ")+1]
                hit_seq_ids[hit_seq_id] = True
                self.log(console, 'HIT: '+hit_seq_id)  # DEBUG
        

        # Filter many set and make filtered output file
        #
        #  Note: don't use SeqIO.parse because loads everything into memory
        #
        self.log(console, 'filtering many sequences')

        if many_type_name == 'SingleEndLibrary':
            
#            with open(many_forward_reads_file_path, 'r', -1) as many_forward_reads_file_handle, open(output_filtered_fasta_file_path, 'w', -1) as output_filtered_fasta_file_handle:
            output_filtered_fasta_file_handle = open(output_filtered_fasta_file_path, 'w', -1)
            if many_forward_reads_file_compression == 'gz':
                many_forward_reads_file_handle = gzip.open(many_forward_reads_file_path, 'r', -1)
            else:
                many_forward_reads_file_handle = open(many_forward_reads_file_path, 'r', -1)

                
            seq_total = 0;
            filtered_seq_total = 0
            last_seq_buf = []
            last_seq_id = None
            last_header = None
            for line in many_forward_reads_file_handle:
                if line.startswith('>'):
                    #self.log(console, 'LINE: '+line)  # DEBUG
                    seq_total += 1
                    seq_id = line[1:]
                    if "\n" in seq_id:
                        seq_id = seq_id[0:seq_id.find("\n")+1]
                    if "\t" in seq_id:
                        seq_id = seq_id[0:seq_id.find("\t")+1]
                    if " " in seq_id:
                        seq_id = seq_id[0:seq_id.find(" ")+1]
                    
                    if last_seq_id != None:
                        #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                        try:
                            in_filtered_set = hit_seq_ids[last_seq_id]
                            #self.log(console, 'HIT')  # DEBUG
                            filtered_seq_total += 1
                            output_filtered_fasta_file_handle.write(last_header)
                            output_filtered_fasta_file_handle.writelines(last_seq_buf)
                        except:
                            pass
                        
                    last_seq_buf = []
                    last_seq_id = seq_id
                    last_header = line
                else:
                    last_seq_buf.append(line)
            if last_seq_id != None:
                #self.log(console, 'ID: '+last_seq_id)  # DEBUG
                try:
                    in_filtered_set = hit_seq_ids[last_seq_id]
                    #self.log(console, 'HIT')  # DEBUG
                    filtered_seq_total += 1
                    output_filtered_fasta_file_handle.write(last_header)
                    output_filtered_fasta_file_handle.writelines(last_seq_buf)
                except:
                    pass
                
                last_seq_buf = []
                last_seq_id = None
                last_header = None

            many_forward_reads_file_handle.close()
            output_filtered_fasta_file_handle.close()

        if filtered_seq_total != hit_total:
            raise ValueError('hits in VSearch alignment output '+str(hit_total)+' != '+str(filtered_seq_total)+' matched sequences in input file')


        #elif many_type_name == 'FeatureSet' \
        #    or many_type_name == 'Genome' \
        #    or many_type_name == 'GenomeSet':
        #
        #    # build a feature set

        # Warning: this reads everything into memory!  Will not work if 
        # the output is very large!
        """
        contigset_data = {
            'id':'megahit.contigset',
            'source':'User assembled contigs from reads in KBase',
            'source_id':'none',
            'md5': 'md5 of what? concat seq? concat md5s?',
            'contigs':[]
        }

        lengths = []
        for seq_record in SeqIO.parse(output_contigs, 'fasta'):
            contig = {
                'id':seq_record.id,
                'name':seq_record.name,
                'description':seq_record.description,
                'length':len(seq_record.seq),
                'sequence':str(seq_record.seq),
                'md5':hashlib.md5(str(seq_record.seq)).hexdigest()
            }
            lengths.append(contig['length'])
            contigset_data['contigs'].append(contig)
        """

        # load the method provenance from the context object
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_one_name'])
        provenance[0]['input_ws_objects'].append(params['workspace_name']+'/'+params['input_many_name'])
        provenance[0]['service'] = 'kb_vsearch'
        provenance[0]['method'] = 'VSearch_BasicSearch'


        # upload reads
        #
        self.upload_SingleEndLibrary_to_shock_and_ws (ctx,
                                                      params['workspace_name'],
                                                      params['output_filtered_name'],
                                                      output_filtered_fasta_file_path,
                                                      provenance
                                                      #sequencing_tech
                                                     )
        # build output report object
        #
        report += 'sequences in many set: '+str(seq_total)
        report += 'sequences in hit set:  '+str(hit_total)

        reportObj = {
            'objects_created':[{'ref':params['workspace_name']+'/'+params['output_filtered_name'], 'description':'SingleEndLibrary VSearch_BasicSearch hits'}],
            'text_message':report
        }

        reportName = 'vsearch_report_'+str(hex(uuid.getnode()))
        report_obj_info = ws.save_objects({
                'id':info[6],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        returnVal = { 'output_report_name': reportName,
                      'output_report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]),
                      'output_filtered_ref': params['workspace_name']+'/'+params['output_filtered_name']
                      }

        #END VSearch_BasicSearch

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method VSearch_BasicSearch return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
