"""
Functions and Classes required for Candida/Mouse intersection notebooks
"""
import os
import re
import requests
import subprocess
import tempfile
import pandas as pd

from Bio import SeqIO
from copy import deepcopy
from yaspin import yaspin

class PepOverlap:
	"Class for storing parsed peptide results"
	def __init__(self) -> None:
		"""
		Creates empty parameters as follows:
			total_dfs(dict) Keyed on species name storing lists of all peptides
			high_freq_dfs(dict) Keyed on species name storing dataframes of peptides present >2 times in the organism
			unique_dfs(list): stores dataframes for each organism of peptides unique within the organism
			summary_data(list) stores per-organism summary dicts
			peptides_total(dict): Keyed on species name; stores lists of all peptides per species
			peptides_unique(dict): Keyed on species name; stores lists of all unique peptides per species
			peptides_high_freq(dict): Keyed on species name; stores lists of all peptides with >2 occurrences in each species
		"""
		self.total_dfs={}
		self.high_freq_dfs={}
		self.unique_dfs=[]
		self.summary_data=[]
		self.peplists_total={}
		self.peplists_unique={}
		self.peplists_high_freq={}
		self.duplicate_df=None

		pass

def download(url,species):
	"""
	Downloads data from specified URL and stores as 'species.fa' in data directory
	Uniprot IDs are reformatted from i.e.

	tr|A0A1D8PCA8|A0A1D8PCA8_CANAL

	to 
	
	A0A1D8PCA8
	
	Required arguments:
		url(str): URL to download
		species(str): species name for output file
	"""
	
	id_re=re.compile('(sp|tr)\|([A-Z0-9\-]*)')
	local_file=f'data/{species}.fa'
	
	if not os.path.exists(local_file):
		tf=next(tempfile._get_candidate_names())
		with yaspin(text=f"Downloading {species} sequences..."):
			with requests.get(url,stream=True) as r:
				r.raise_for_status()
				with open(f'data/{tf}','wb') as fh:
					for dat in r.iter_content(8192):
						fh.write(dat)
	
		#Rewrite sequence IDs in simpler format for EMBOSS, and clean up description
		with open(local_file,'w') as out_fh:
			with open(f'data/{tf}','r') as in_fh:
				for record in SeqIO.parse(in_fh,'fasta'):
					match=id_re.search(record.id)
					if match:
						record.id=match.group(2)
					else: 
						raise Exception(f'Failed to parse identifier: {record.id}')
					
					desc_fields=record.description.split(' ')
					desc_fields.pop(0)
					record.description=' '.join(desc_fields)
					
					SeqIO.write([record],out_fh,'fasta')
						
		os.remove(f'data/{tf}')

def pepdigest(beast):
	"""
	Runs EMBOSS pepdigest on specified species, assuming fasta file 'species.fa' available in data directory
	
	Required params:
		beast(str): organism name
	
	Returns:
		None
	"""
	if not os.path.exists(f'data/{beast}.pepdigest'):
		with yaspin(text=f"Running {beast} pepdigest"):
			subprocess.run(["pepdigest", "-auto", "-seqall", f'data/{beast}.fa',
							'-menu', '1', '-mono','N','-outfile',f'data/{beast}.pepdigest'],
						   check=True)

	
def parse_pepdigest(beast): 
	
	"""
	Parses pepdigest output file storing resulting data in pandas dataframe
	
	Required params:
		beast(str): organism name
		
	Returns:
		df(pd.DataFrame): Pandas dataframe with one row per peptide. The following columns
						  are present in the dataframe:
						  
						  protein_id: Uniprot protein accession
						  count:      No. of occurrences of peptide in proteome
						  peptide:    AA sequence of peptide
						  frequency:  No. of occurrences of peptide in proteome
	"""
	
	digestfile=f'data/{beast}.pepdigest'

	# new_seq_re matches the start of each record beginning capturing the uniprot ID in group 1 i.e.
	# # Sequence: A0A087WPF7     from: 1   to: 1261
	new_seq_re = re.compile('# Sequence: ([\S]+)')

	# peptide_re matches each line containing a peptide, capturing the peptide sequence in group 3 i.e.
	# Start	End	Mol_Weight	Cterm	Nterm	Sequence
	# 57	78	2541.862	K		E		EDNGKPPSSAPSRPRPPRRKRR
	peptide_re = re.compile('([\s]+[\d\.]+){3}[\s]+([A-Z\.][\s]+){2}([A-Z]+)')
	
	with open(digestfile, 'r') as fh:
		prot_id=None
		pep_count=0
		peptides=[]
	
		for line in fh:
			match = new_seq_re.match(line)
			if match:
				prot_id=match.group(1)
				pep_count=0
			
			match = peptide_re.match(line)
			if match:
				pep_count=pep_count+1
				peptide=match.group(3) 
				pep_data={
					'protein_id': prot_id,
					'count': pep_count,
					'peptide': peptide,
					'organism': beast
				}
				peptides.append(pep_data)
		
	columns=['protein_id','count','peptide']
	df = pd.DataFrame(peptides,columns=columns)

	# Some peptides are present multipe times, so determine the frequency of each peptide 
	# occurring and add to the dataframe
	df['frequency']=df.groupby('peptide')['peptide'].transform('count')

	return df

def summarise_proteome(beast,pep_data):
	"""
	Produces summary data for a pepdigest dataframe
	
		df(pd.DataFrame): organism pepdigest parsed data as returned by parse_pepdigest()
		beast(str): organism name
		pep_data(PepOverlap): Object for storing parsed peptide info
		
	returns:
		pep_data(PepOverlap)
	"""
	
	df=pep_data.total_dfs[beast]
	
	unique_df=df.drop_duplicates('peptide',keep=False)
	unique_df=deepcopy(unique_df)
	unique_df.loc[:,'organism']=beast
	pep_data.unique_dfs.append(unique_df)
	
	# high_freq refers to peptides present more than two times in the proteome
	high_freq_df=df[df['frequency']>2]
	high_freq_df=high_freq_df.sort_values(by=['frequency'])
	pep_data.high_freq_dfs[beast]=high_freq_df
	
	# peplist dicts store a list of each class of peptide (total, unique, high_freq)
	pep_data.peplists_total[beast]=df['peptide'].tolist()
	pep_data.peplists_unique[beast]=unique_df['peptide'].tolist()
	pep_data.peplists_high_freq[beast]=high_freq_df['peptide'].tolist()

	# Now collect some summary statistics, to be stored in the 'summary_data' list
	protein_count=len(df['protein_id'].unique().tolist())
	unique_protein_count=len(unique_df['protein_id'].unique().tolist())
	peptide_count=(len(df.index))
	unique_peptide_count=(len(unique_df.index))

	data={
		'Category': beast,
		'Peptide count': peptide_count,
		'Protein count': protein_count,
		'Unique peptides': unique_peptide_count,
		'Unique proteins': unique_protein_count
	}
	pep_data.summary_data.append(data)

	return(pep_data)
	
def get_prot_info(proteins):
	"""
	Parses metadata from SeqRecords (assuming Uniprot fasta-format description lines)
	
	Required parameters:
		proteins(list) - list of Bio SeqRecords
		
	Returns:
		df(pd.DataFrame) - dataframe with 1 row per protein.  The following columns
						   are present in the dataframe:
						  
						   Accession: Uniprot protein accession
						   Description: Description from Uniprot record
						   Gene Symbol: Gene Symbol from Uniprot record
	"""
	info=[]
	symbol_re=re.compile('GN=([\S]+)')
	desc_re=re.compile('[A-Z0-9-]* ([^$]+) OS=')
	
	for prot in proteins:
		match=symbol_re.search(prot.description)
		# Not all entries include a gene symbol...
		if match:
			gene_symbol=match.group(1)
		else:
			gene_symbol=''
		
		match=desc_re.search(prot.description)
	
		if match:
			description=match.group(1)
		else:
			raise Exception('no description in {}'.format(desc))
	
		info.append({'Accession'  : prot.id,
					 'Description': description,
					 'Gene Symbol': gene_symbol})
	
	columns=('ID', 'Gene Symbol','Description')
	df=pd.DataFrame(info)

	return(df)

def format_venn(plot):
	"""
	Apply standardised formatting to venn diagrams
	
	Required parameters:
		plot(venn2): venn2 plot

	Returns:
		plot(venn2): reformatted plot
	"""
	for text in plot.set_labels:
		text.set_fontsize(14)
	for text in plot.subset_labels:
			text.set_fontsize(14)
	plot.get_patch_by_id('01').set_color('#62b3e6')
	plot.get_patch_by_id('01').set_edgecolor('none')
	plot.get_patch_by_id('11').set_color('#464646')
	plot.get_patch_by_id('11').set_edgecolor('none')
	plot.get_patch_by_id('10').set_color('#62b3e6')
	plot.get_patch_by_id('10').set_edgecolor('none')
	if plot.get_label_by_id('01').get_text()=='0':
		plot.get_label_by_id('01').set_text('')

	return(plot)

def clean_peptide(peptide):
	"""
	Clean the peptide sequence from the spectronaut outputs to remove anything
	other than peptide sequence...

	Required params:
	  peptide(str)

	Returns:
	  peptide(str)
	"""
	# strip out flanking '_' and any trailing '.[0-9]'
	peptide=re.sub('_(\.[0-9]+)?','',peptide)
	# remove any modification info which is contained within square brackets
	peptide=re.sub('\[[A-Za-z\(\)\- ]*\]','',peptide)

	return(peptide)

def parse_fasta_desc(filename):
	"""
	Parses metadata from Uniprot fasta formatted description lines

	Requred parameters:
		filename(str): path to fasta file
	
	Returns:
		df(pd.DataFrame): Pandas DataFrame 
	"""
	desc_re=re.compile('[A-Z0-9-]* (.*) OS=[A-Za-z0-9\/\-\(\) ]*OX=[0-9]*( GN=([A-Za-z0-9]*))?')

	descs=list()
	with open(filename,'r') as fh:
		for record in SeqIO.parse(fh,'fasta'):
			match=re.match(desc_re,record.description)
			if match:
				descs.append({
					'prot_id': record.id,
					'description': match.group(1),
					'gene': match.group(3)
				})
			else:
				raise Exception(f'Failed to parse description: {record.description}')

	df=pd.DataFrame(descs)

	return(df)